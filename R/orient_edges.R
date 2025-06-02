
# function: for all nodes, fit GAM using parents and returns residuals from fit
# input: data, adjacency matrix, name of nodes, number of basis functions (k), whether to count neighboring nodes as parents
# output: matrix containing the residuals
gam.residuals <- function(data, adjmat, nodelabels, k, use_nbrs=FALSE){
  dat_rsd = data

  # gam parameters
  gam.control.param = list(nthreads=4, epsilon=1e-07, maxit=150, mgcv.tol=1e-07)

  for(i in 1:ncol(data)) {
    gc()
    # parents of node i
    if(use_nbrs==TRUE){
      pa_nodes = which(adjmat[, i] == 1)
    }else{
      pa_nodes = which(adjmat[, i] - adjmat[i, ] == 1)
    }
    # residuals of node i
    if(sum(pa_nodes)) {
      fm = paste(nodelabels[i], "~", paste("s(", pa_nodes, ", k=", k, ", bs='tp')",
                                           sep = "", collapse = " + "))
      gfit = mgcv::gam(formula = as.formula(fm), data = as.data.frame(data), optimizer = c("outer","newton"),
                 select = TRUE, gam.control=gam.control.param, debug=TRUE)
      dat_rsd[, i] = gfit$residuals
    }
  }
  return(dat_rsd)
}

# function: fit GAM
# input: data, covariates, response variable, number of basis functions (k)
# output: model object
fit_gam_reg <- function(data, pa_nodes, ch_node, k){
  gam.control.param = list(nthreads=4, epsilon=1e-07, maxit=150, mgcv.tol=1e-07)
  if(length(pa_nodes) > 0) {
    fm = paste(ch_node, "~", paste("s(", pa_nodes, ", k=", k, ", bs='tp')",
                                   sep = "", collapse = " + "))
    main_gfit = mgcv::gam(formula = as.formula(fm), data = as.data.frame(data), optimizer = c("outer","newton"),
                    select = TRUE, gam.control=gam.control.param, debug=TRUE)
  }else{
    # main_gfit <- data[,ch_node]
    main_gfit <- lm(data[,ch_node]~1)
  }
  return(main_gfit)
}

# function: fit GAM each a pair of nodes to produce PANM
# input: data, adjacency matrix, name of nodes, covariates, response variable, number of basis functions (k)
# output: 2 model objects in list
gam.pa.fit <- function(data, adjmat, nodelabels, pa_node, ch_node, k){
  # parents of child node, including candidate parent
  if(is.numeric(pa_node)){
    pa_node <- nodelabels[pa_node]
  }
  pa_nodes = c(nodelabels[unname(which(adjmat[, ch_node] - adjmat[ch_node, ] == 1))], pa_node)

  if(is.numeric(ch_node)){
    ch_node <- nodelabels[ch_node]
  }

  main_gfit <- fit_gam_reg(data, pa_nodes, ch_node, k)

  # parents of candidate parent node
  pa_nodes_prior = nodelabels[which(adjmat[, pa_node] - adjmat[pa_node, ] == 1)]
  prior_gfit <- fit_gam_reg(data, pa_nodes_prior, pa_node, k)
  return(list(main_gfit, prior_gfit))
}


# function: check if edge orientation induces cycle in PDAG
# input: adjacency matrix, child node, parent node, whether to assess as CPDAG or DAG
# output: true/false

#' Check if adding edge induces cycle in graph
#'
#' @param adjmat adjacency matrix of graph
#' @param ch_ind name/index of child node
#' @param pa_ind name/index of parent node
#' @param cpdag evaluate as CPDAG or not (if not, evaluated as DAG)
#'
#' @return TRUE if edge induces cycle, FALSE otherwise
#' @export
#'
#' @examples
#' creates_cycles(amat, "X1", "X2", cpdag=T)
creates_cycles <- function(adjmat, ch_ind, pa_ind, cpdag=T){
  # a node is needed.
  if(missing(adjmat)){
    stop("no adjacency matrix specified.")
  }

  adjmat[ch_ind, pa_ind] <- 0
  adjmat[pa_ind, ch_ind] <- 1

  temp_dag <- bnlearn::empty.graph(nodes=colnames(adjmat))
  amat(temp_dag, check.cycles=FALSE) <- adjmat


  # check whether the PDAG contains completely directed cycles
  # CPDAG represents whether we are assessing a CPDAG or DAG
  if(bnlearn::acyclic(temp_dag, directed = cpdag)){
    return(FALSE) #return true if it contains cycles
  }else{
    return(TRUE)
  }

}


# function: count number of neighbors(connected by undirected edge) shared by two nodes
# input: adjacency matrix, vector containing names of nodes
# output: number of shared neighbors

count_common_nbr_ind <- function(adjmat, curr_edge){
  n1 <- curr_edge[1]
  n2 <- curr_edge[2]

  # get undirected edges for each node
  n1_nbr <- adjmat[n1, ] == adjmat[, n1] & adjmat[n1, ] == 1
  n2_nbr <- adjmat[n2, ] == adjmat[, n2] & adjmat[n2, ] == 1
  # count common neighbor
  return(sum(n1_nbr + n2_nbr == 2))
}

# function: count number of neighbors(connected by undirected edge) shared for all undirected edges
# input: adjacency matrix, matrix of undirected edges
# output: vector containing number of shared neighbors for each edge

count_common_nbr <- function(adjmat, udr_edges, existing_nbr=NULL){
  if(is.null(nrow(udr_edges))){
    return(c(0))
  }
  n_common_nbr <- rep(0, nrow(udr_edges))
  # if the number of common nbr is unknown for all edges
  # this is meant for getting initial counts
  if(is.null(existing_nbr)){
    for(i in c(1:nrow(udr_edges))){
      n_common_nbr[i] <- count_common_nbr_ind(adjmat, curr_edge=udr_edges[i,])
    }
  }else{
    # if number of common nbr known for some
    # this is meant for updating the counts
    for(i in c(1:nrow(udr_edges))){
      if(existing_nbr[i] > 0){
        n_common_nbr[i] <- count_common_nbr_ind(adjmat, curr_edge=udr_edges[i,])
      }
    }
  }
  return(n_common_nbr)
}


# function: compute the standardized mutual information between 2 variables
# input: data of var1, data of var2, discretization method, number of bins
# output: vector containing number of shared neighbors for each edge
#' Compute the standardized mutual information between 2 variables
#'
#' @param var1 vector containing data of variable 1
#' @param var2 vector containing data of variable 2
#' @param disc_method discretization method
#' @param nbins number of bins
#'
#' @return (numeric value) standardized mutual information
#' @export
#'
#' @examples
#' compute_mi_std(v1, v2)
compute_mi_std <- function(var1, var2, disc_method='equalwidth', nbins=20){
  disc_var1 <- as.vector(infotheo::discretize(var1,disc=disc_method, nbins=nbins)$X)
  disc_var2 <- as.vector(infotheo::discretize(var2,disc=disc_method, nbins=nbins)$X)
  var1_entropy <- infotheo::entropy(disc_var1)
  var2_entropy <- infotheo::entropy(disc_var1)
  mi_std <- infotheo::mutinformation(X=disc_var1, Y=disc_var2)/min(var1_entropy, var2_entropy)
  return(mi_std)
}


# function: computes the independence measure (adherence to PANM assumption) for an edge
# input: data, adj matrix, vector containing nodes of edge, number of basis functions in GAM (k), train/test split ratio
# output: std. MI for X->Y, std. MI for Y->X

compute_panm_indep <- function(data, adjmat, curr_edge, k, train_ratio=0.5,
                               compute_min_indep_only=FALSE, train_index=NULL, samples_per_bin=50){
  nodelabels <- colnames(data)
  n1 <- curr_edge[1]
  n2 <- curr_edge[2]
  if(is.numeric(n1)) n1 <- nodelabels[n1]
  if(is.numeric(n2)) n2 <- nodelabels[n2]
  # train/test split on data
  if(is.null(train_index)){
    train_index <- sample(c(1:nrow(data)), size=nrow(data)*train_ratio, replace = FALSE)
  }
  train_data <- as.data.frame(base::apply(data[train_index,], MARGIN = 2, FUN = base::scale))
  test_data <- as.data.frame(base::apply(data[-train_index,], MARGIN = 2, FUN = base::scale))

  # test independence between residual of n1 & n2, then residual of n2 & n1
  # if both results show dependence, then n1 & n2 have a common parent
  # build regression model for n1 = pa(n1) + n2
  n1_pa <- nodelabels[which(adjmat[, n1] - adjmat[n1, ] == 1)]
  model_n1 <- fit_gam_reg(data=train_data, pa_nodes=c(n1_pa, n2), ch_node=n1, k=k)
  # predict n1 on testing data and get residuals
  n1_hat <- predict(model_n1, newdata=as.data.frame(test_data), type='response')
  n1_predictor <- predict(model_n1, newdata=as.data.frame(test_data), type='iterms')[,paste("s(", n2, ")", sep="")]
  n1_residuals <- test_data[,n1] - n1_hat

  # check other node
  n2_pa <- nodelabels[which(adjmat[, n2] - adjmat[n2, ] == 1)]
  model_n2 <- fit_gam_reg(data=train_data, pa_nodes=c(n2_pa, n1), ch_node=n2, k=k)
  # predict n1 on testing data and get residuals
  n2_predictor <- predict(model_n2, newdata=as.data.frame(test_data), type='iterms')[,paste("s(", n1, ")", sep="")]
  n2_hat <- predict(model_n2, newdata=as.data.frame(test_data), type='response')
  n2_residuals <- test_data[,n2] - n2_hat

  # discretize data
  disc_method <- "equalwidth"
  nbins <- round(nrow(test_data)/samples_per_bin)

  # compute standardized mi for res(n1) and n2
  mi_n1_std <- compute_mi_std(n1_residuals, test_data[,n2], disc_method=disc_method, nbins=nbins)
  # compute standardized mi for res(n2) and n1
  mi_n2_std <- compute_mi_std(n2_residuals, test_data[,n1], disc_method=disc_method, nbins=nbins)
  # compute standardized mi for res(n1) and parents of n1
  mi_n1_std_pa <- c()
  if(length(n1_pa)==0 | is.null(n1_pa)){
    mi_n1_std_pa <- mi_n1_std
  }else{
    for(pa_node in n1_pa){
      mi_pa <- compute_mi_std(n1_residuals, test_data[,pa_node])
      mi_n1_std_pa <- c(mi_n1_std_pa, mi_pa)
    }
  }
  mi_n1_std_pa <- max(mi_n1_std_pa)

  # compute standardized mi for res(n2) and parents of n2
  mi_n2_std_pa <- c()
  if(length(n2_pa)==0 | is.null(n2_pa)){
    mi_n2_std_pa <- mi_n2_std
  }else{
    for(pa_node in n2_pa){
      mi_pa <- compute_mi_std(n2_residuals, test_data[,pa_node])
      mi_n2_std_pa <- c(mi_n2_std_pa, mi_pa)
    }
  }
  mi_n2_std_pa <- max(mi_n2_std_pa)

  if(compute_min_indep_only == TRUE){
    return(min(mi_n1_std, mi_n2_std, mi_n1_std_pa, mi_n2_std_pa))
  }

  return(c(mi_n1_std, mi_n2_std))
}


#' Perform likelihood ratio based orientation test on an undirected edge
#'
#' Corresponds to the pairwise Additive Noise Models built for each direction in the paper
#'
#' @param data observed data
#' @param train_ratio ratio ranging from 0-1 for train/test split on data
#' @param ll_test_approach "2fold" for running the 2 fold test, other for 5-5 train test split approach
#' @param adjmat adjacency matrix (of maximally oriented CPDAG)
#' @param pa_ind index/name of parent node
#' @param ch_ind index/name of child node
#' @param k_basis number of basis functions for GAM
#'
#' @return list containing likelihood ratio test result and statistics
#' @examples
#' ll_orient_test(data, ll_test_approach='55', adjmat=adjmat, pa_ind=2, ch_ind=5)
#' @export

ll_orient_test <- function(data, train_ratio=0.5, ll_test_approach="SS", adjmat, pa_ind, ch_ind, k_basis){
  # train/test split
  train_index <- sample(c(1:nrow(data)), size=nrow(data)*train_ratio, replace = FALSE)
  train_data <- as.data.frame(base::apply(data[train_index,], MARGIN = 2, FUN = base::scale))
  test_data <- as.data.frame(base::apply(data[-train_index,], MARGIN = 2, FUN = base::scale))

  if(ll_test_approach=="CV"){
    # parent node is the node in col 1
    # first fold
    gam_XtoY_m1 <- gam.pa.fit(train_data, adjmat, nodelabels=colnames(data), pa_node=pa_ind, ch_node=ch_ind, k=k_basis)
    gam_YtoX_m1 <- gam.pa.fit(train_data, adjmat, nodelabels=colnames(data), pa_node=ch_ind, ch_node=pa_ind, k=k_basis)
    # second fold
    gam_XtoY_m2 <- gam.pa.fit(test_data, adjmat, nodelabels=colnames(data), pa_node=pa_ind, ch_node=ch_ind, k=k_basis)
    gam_YtoX_m2 <- gam.pa.fit(test_data, adjmat, nodelabels=colnames(data), pa_node=ch_ind, ch_node=pa_ind, k=k_basis)
    arc_est <- get_edge_dir_2fold(gam_XtoY_m1, gam_YtoX_m1, gam_XtoY_m2, gam_YtoX_m2,
                                  train_data1=train_data, train_data2=test_data)
  }else{
    # 5/5 split, one test
    gam_XtoY <- gam.pa.fit(train_data, adjmat, nodelabels=colnames(data), pa_node=pa_ind, ch_node=ch_ind, k=k_basis)
    gam_YtoX <- gam.pa.fit(train_data, adjmat, nodelabels=colnames(data), pa_node=ch_ind, ch_node=pa_ind, k=k_basis)
    arc_est <- get_edge_dir(gam_XtoY[[1]], gam_YtoX[[1]], gam_XtoY[[2]], gam_YtoX[[2]], test_data = test_data)
  }
  return(arc_est)
}


# function: orient edge in the PDAG
# input: pdag (bnlearn object), adjacency matrix, index of parent node, index of child node
# output: pdag (bnlearn object) with edge oriented and added to whitelist
update_bnlearn_dag_wl <- function(bnlearn_obj, pa_ind, ch_ind){
  adjmat <- amat(bnlearn_obj)
  nodelabels <- colnames(adjmat)
  if(!is.null(adjmat)){
    if(adjmat[ch_ind, pa_ind]==1){
      temp_ind <- ch_ind
      ch_ind <- pa_ind
      pa_ind <- temp_ind
    }
  }
  if(is.null(bnlearn_obj$learning$whitelist)){
    wl <- matrix(c(nodelabels[pa_ind], nodelabels[ch_ind]), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
    bnlearn_obj$learning$whitelist <- wl
  }else{
    bnlearn_obj$learning$whitelist <- rbind(bnlearn_obj$learning$whitelist, c(nodelabels[pa_ind], nodelabels[ch_ind]))
  }
  return(bnlearn_obj)
}


#' STEP 2 OF ALGORITHM: orient undirected edges in the CPDAG
#'
#' The procedure first categorizes edges by number of neighbors shared between nodes,
#' Computes an independence measure for each edge, which is used to sort edges,
#' And performs the likelihood test to determine the causal direction
#'
#' @param data observed data
#' @param curr_amat adjacency matrix of graph object
#' @param alpha significance level for likelihood test
#' @param ll_test_approach "2fold" for running the 2 fold test, other for 5-5 train test split approach
#' @param nodelabels vector of variable names
#' @param train_ratio percentage of data used for training data
#'
#' @return adjacency matrix of updated PDAG
#' @examples
#' orient_edge(data, curr_amat=adjmat, ll_test_approach='55')
#' @export
orient_edge = function(data, curr_amat, nodelabels = colnames(data), alpha = 0.05,
                       train_ratio=0.5, ll_test_approach=NULL){
  # extract adjacency matrix
  adjmat = curr_amat
  curr_cpdag = bnlearn::empty.graph(colnames(curr_amat))
  amat(curr_cpdag) = curr_amat

  # check if this is already a DAG; screen for undirected edges
  if(sum(adjmat + t(adjmat) == 2) == 0){
    out = list("amat"= adjmat, "cpdag" = curr_cpdag, "sigpair" = NULL,
               "full_res" = NULL, "total_time" = 0, "edgeorient_time" = 0)
    return(out)
  }

  if(is.null(nodelabels)){
    nodelabels <- paste("X", 1:ncol(data), sep = "")
    colnames(data) = nodelabels
    colnames(adjmat) = nodelabels
  }
  # assign test method based on sample size
  if(is.null(ll_test_approach)) ll_test_approach <- ifelse(nrow(data)>500, "2fold", "5_5_split")

  total_start = Sys.time()
  # undirected arcs and corresponding pairs of nodes
  udr_edges = matrix(which((adjmat == t(adjmat) & (adjmat != 0)), arr.ind = TRUE), ncol = 2)
  # just keep one record of the edges; previous line returns both directions ie two rows
  udr_edges = matrix(udr_edges[udr_edges[, 1] < udr_edges[, 2], ], ncol = 2)

  # calculate and sort by number of common neighbors
  udr_edges <- cbind(udr_edges, count_common_nbr(adjmat, udr_edges=udr_edges))
  udr_edges <- udr_edges[order(udr_edges[,3], decreasing=FALSE),,drop=FALSE]
  max_common_nbr <- max(udr_edges[,3])

  # adding column for storing p-values, diff in ll, and whether edge is oriented
  # using 1 & 0 as placeholders
  udr_edges <- cbind(udr_edges, rep(1, nrow(udr_edges)), rep(0, nrow(udr_edges)), rep(0, nrow(udr_edges)), rep(0, nrow(udr_edges)))
  colnames(udr_edges) <- c("pa", "ch", "n_nbr", "pval", "ll_diff", "has_eval", "has_orient")

  # for saving results
  res = vector("list")
  sig_ind = NULL
  ll_start <- Sys.time()

  # number of base functions for smoothing function
  if(nrow(data)/ncol(data) < 3*10) k_basis <- ceiling(nrow(data)/(3*ncol(data)))
  else k_basis <- min(7, ceiling(sqrt(ncol(data))-1))

  ll_test_start <- Sys.time()

  # start with neighors n_nbrs=0,1,2,..
  # for edges with common nbr(s), rank by mutual information
  for(curr_n_nbr in c(0:max_common_nbr)){
    udr_edges_subset <- udr_edges[udr_edges[,"n_nbr"]==curr_n_nbr,,drop=F]
    if(nrow(udr_edges_subset)==0) next
    # if there are common nbrs, rank edges
    # if(curr_n_nbr > 0 && nrow(udr_edges_subset > 1)){
    if(nrow(udr_edges_subset > 1)){
      min_mi_vec <- apply(udr_edges_subset, MARGIN=1,
                          FUN=function(x) has_common_pa(dat, adjmat=adjmat, curr_edge=x,
                                                        k=k_basis, compute_min_indep_only = FALSE))
      new_colnames <- c(colnames(udr_edges_subset), "mi_n1_std_pa", "mi_n2_std_pa", "min_mi")
      udr_edges_subset <- cbind(udr_edges_subset, t(min_mi_vec))
      # colnames(udr_edges_subset)[ncol(udr_edges_subset)] <- "min_mi"
      colnames(udr_edges_subset) <- new_colnames
      # udr_edges_subset <- udr_edges_subset[order(udr_edges_subset[, "min_mi"], decreasing=FALSE),,drop=FALSE]
    }


    for(i in c(1:nrow(udr_edges_subset))){
      # get parent and child indices, with row num in larger list
      pa_ind = udr_edges_subset[i, "pa"]
      ch_ind = udr_edges_subset[i, "ch"]
      orig_df_index <- which(udr_edges[,"pa"] == pa_ind & udr_edges[,"ch"] == ch_ind, arr.ind=TRUE)
      udr_edges[orig_df_index, "has_eval"] <- 1

      # CHECK IF EDGE BECOMES DIRECTED AFTER EDGE ORIENTATION, SKIP TO SAVE TIME
      if(adjmat[nodelabels[pa_ind], nodelabels[ch_ind]]-adjmat[nodelabels[ch_ind],nodelabels[pa_ind]]!=0){
        udr_edges[orig_df_index, "has_orient"] <- 1
        gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, ch_ind)
        next
      }

      # perform likelihood test
      ll_test_res <- ll_orient_test(dat, train_ratio=0.5, ll_test_approach=ll_test_approach,
                                    adjmat=adjmat, pa_ind=pa_ind, ch_ind=ch_ind, k_basis=k_basis)
      res <- append(res, ll_test_res)
      udr_edges[orig_df_index, "pval"] <- ll_test_res$pval_ll_diff
      udr_edges[orig_df_index, "ll_diff"] <- ll_test_res$llratio

      # causal dirction found to be significant
      if(ll_test_res$pval_ll_diff < alpha){
        # change parent node if second model is better
        if(ll_test_res$llratio < 0){
          temp_node <- pa_ind
          pa_ind <- ch_ind
          ch_ind <- temp_node
        }
        # check for cycles
        if(creates_cycles(adjmat, ch_ind, pa_ind)==TRUE) next

        # update existing CPDAG and orient edges again
        gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[pa_ind], to=nodelabels[ch_ind])
        adjmat = amat(gam_cpdag)
        sig_ind = rbind(sig_ind, c(pa_ind, ch_ind))
        udr_edges[orig_df_index, "has_orient"] <- 1
        # update whitelist in cpdag
        gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, ch_ind)

        # orient all edges of common neighbors & children outwards
        pa_node_nbr <- nodelabels[adjmat[,nodelabels[pa_ind]] == adjmat[nodelabels[pa_ind],] & adjmat[,nodelabels[pa_ind]] == 1]
        pa_node_ch <- nodelabels[adjmat[nodelabels[pa_ind],] - adjmat[,nodelabels[pa_ind]] == 1]
        ch_node_nbr <- nodelabels[adjmat[,nodelabels[ch_ind]] == adjmat[nodelabels[ch_ind],] & adjmat[,nodelabels[ch_ind]] == 1]
        ch_node_ch <- nodelabels[adjmat[nodelabels[ch_ind],] - adjmat[,nodelabels[ch_ind]] == 1]

        union_pa_nc <- union(pa_node_nbr, pa_node_ch)
        union_ch_nc <- union(ch_node_nbr, ch_node_ch)

        common_nc <- intersect(union_pa_nc, union_ch_nc)
        common_nc = common_nc[!(common_nc %in% nodelabels[c(pa_ind, ch_ind)])]

        # orient edges as inward for all common neighbors & children
        if(length(common_nc) > 0){
          for(nbr_ind in common_nc){
            if(creates_cycles(adjmat, nbr_ind, pa_ind)==TRUE || creates_cycles(adjmat, nbr_ind, ch_ind)==TRUE) next
            gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[pa_ind], to=nbr_ind)
            gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[ch_ind], to=nbr_ind)
            adjmat = amat(gam_cpdag)
            # update whitelist in cpdag
            gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, which(nodelabels == nbr_ind))
            gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, ch_ind, which(nodelabels == nbr_ind))
          }
        }

        # apply Meek's rule (ie convert to cpdag) to orient more edges
        gam_cpdag <- cpdag(gam_cpdag, wlbl=TRUE)
        adjmat = amat(gam_cpdag)

        # # update whitelist in cpdag
        # gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, adjmat=NULL, pa_ind, ch_ind)

        # update number of neighbors
        udr_edges[,"n_nbr"] <- count_common_nbr(adjmat, udr_edges[,c(1,2)], existing_nbr=udr_edges[,"n_nbr"])
      }else{next}
      gc()
    }
  }

  # orient remaining edges
  udr_edges_subset <- NULL
  if(sum(udr_edges[,"has_eval"] == 0) + sum(udr_edges[,"pval"]<alpha & udr_edges[,"has_orient"]==0) > 0){
    # edges not evaluated + edge found sig but weren't oriented
    udr_edges_subset <- udr_edges[(udr_edges[,"has_eval"]==0) | (udr_edges[,"pval"]<alpha & udr_edges[,"has_orient"]==0),,drop=FALSE]
    udr_edges_subset <- udr_edges_subset[order(udr_edges_subset[,"n_nbr"], decreasing=FALSE),,drop=FALSE]
    for(i in c(1:nrow(udr_edges_subset))){
      # get parent and child indices, with row num in larger list
      pa_ind = udr_edges_subset[i, "pa"]
      ch_ind = udr_edges_subset[i, "ch"]
      orig_df_index <- which(udr_edges[,"pa"] == pa_ind & udr_edges[,"ch"] == ch_ind, arr.ind=TRUE)
      udr_edges[orig_df_index, "has_eval"] <- 1

      # CHECK IF EDGE BECOMES DIRECTED AFTER EDGE ORIENTATION, SKIP TO SAVE TIME
      if(adjmat[nodelabels[pa_ind], nodelabels[ch_ind]]-adjmat[nodelabels[ch_ind],nodelabels[pa_ind]]!=0){
        udr_edges[orig_df_index, "has_orient"] <- 1
        gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, ch_ind)
        next
      }

      # perform likelihood test
      ll_test_res <- ll_orient_test(dat, train_ratio=0.5, ll_test_approach=ll_test_approach,
                                    adjmat=adjmat, pa_ind=pa_ind, ch_ind=ch_ind, k_basis=k_basis)
      res <- append(res, ll_test_res)
      udr_edges[orig_df_index, "pval"] <- ll_test_res$pval_ll_diff
      udr_edges[orig_df_index, "ll_diff"] <- ll_test_res$llratio

      # causal dirction found to be significant
      if(ll_test_res$pval_ll_diff < alpha){
        # change parent node if second model is better
        if(ll_test_res$llratio < 0){
          temp_node <- pa_ind
          pa_ind <- ch_ind
          ch_ind <- temp_node
        }
        # check for cycles
        if(creates_cycles(adjmat, ch_ind, pa_ind)==TRUE) next

        # update existing CPDAG and orient edges again
        gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[pa_ind], to=nodelabels[ch_ind])
        sig_ind = rbind(sig_ind, c(pa_ind, ch_ind))
        udr_edges[orig_df_index, "has_orient"] <- 1

        # update whitelist in cpdag
        gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, ch_ind)
        # apply Meek's rule (ie convert to cpdag) to orient more edges
        gam_cpdag <- cpdag(gam_cpdag, wlbl=TRUE)
        adjmat = amat(gam_cpdag)
        # # update number of neighbors
        # udr_edges[,"n_nbr"] <- count_common_nbr(adjmat, udr_edges[,c(1,2)], existing_nbr=udr_edges[,"n_nbr"])
      }else{next}
      gc()
    }
  }

  ll_test_end <- Sys.time()
  if(!is.null(sig_ind)){
    colnames(sig_ind) <- c("from", "to")
    sig_ind[,1] <- nodelabels[as.numeric(sig_ind[,1])]
    sig_ind[,2] <- nodelabels[as.numeric(sig_ind[,2])]
  }

  ll_end <- Sys.time()
  total_time <- as.numeric(ll_end-total_start, unit = "secs")
  edgeorient_time <- as.numeric(ll_test_end-ll_test_start, unit = "secs")
  return(amat(curr_cpdag))
}

