

# purpose: gets parents (directed and undirected edges) of node
# input: skeleton from pPC_skeleton()/list of nodes from bnlearn obj, node
# output: directed and undirected parents of node (with labels of being pa or not) or FALSE if no parents/node
get_curr_pa = function(x, node) {
  # check if node is in skeleton
  if(node %in% names(x)){
    pa_set <- x[[node]]$nbr[!x[[node]]$nbr %in% x[[node]]$children]
    is_confirmed_pa <- pa_set %in% x[[node]]$parents
  }else{
    return(NULL)
  }

  if(length(pa_set) > 0){
    return(cbind(pa_set, is_confirmed_pa))
  }
  else{
    return(NULL)
  }
}

# purpose: removes edge from skeleton & updates properties
# input: PC skeleton, nodes to update, p-value of coef from GAM
# output: skeleton (modifies skeletons then returns it)
drop_arc = function(dag, node1, node2, p_value, dsep_set){
  # update neighbor & markov blanket for node1, node2
  dag[[node1]]$nbr <- dag[[node1]]$nbr[dag[[node1]]$nbr != node2]
  dag[[node2]]$nbr <- dag[[node2]]$nbr[dag[[node2]]$nbr != node1]
  dag[[node1]]$mb <- dag[[node1]]$mb[dag[[node1]]$mb != node2]
  dag[[node2]]$mb <- dag[[node2]]$mb[dag[[node2]]$mb != node1]

  # update p-value & dsep.set in attributes
  if(is.null(attr(dag, "dsep.set"))){
    return(dag)
  }

  p_value <- unname(p_value)
  node_pairs <- attr(dag, "dsep.set")
  pair_to_drop <- sort(c(node1, node2))
  for(i in 1:length(node_pairs)){
    if(all(pair_to_drop == sort(node_pairs[[i]]$arc))){
      attr(dag, "dsep.set")[[i]]$p.value <- p_value
      # add found conditional indep. relation to current dsep_set
      attr(dag, "dsep.set")[[i]]$dsep.set <- c(attr(dag, "dsep.set")[[i]]$dsep.set, dsep_set)
    }
  }
  return(dag)
}


#' Builds GAM formula to be used in gam()
#'
#' @param y response variable name
#' @param x vector containing all covariate names
#' @param k number of basis functions for smoothing splines in GAM
#' @param smooth_fx basis function used
#'
#' @return
#' @export
#'
#' @examples
#' get_gam_formula("Y", c("X", "W", "Z"), k=3, smooth_fx="cr")
get_gam_formula <- function(y, x, k, smooth_fx="tp"){
  covar <- "~"
  for(i in 1:length(x)){
    if(length(x)==0){
      break
    }
    covar <- paste(covar, "s(", x[i], ", k=", k, ", bs='", smooth_fx, "') +", sep="")
  }
  covar <- substr(covar, 1, nchar(covar)-1)
  fm <- as.formula(paste(y, covar, sep=""))
  return(fm)
}



#' STEP 3 OF ALGORITHM: Delete superfluous edges from PDAG to result in DAG
#'
#' For each node, fit regression model using parents and neighboring nodes
#' Perform significance testing on covariates to remove edges
#'
#' @param data observed data
#' @param curr_amat adjacency matrix of graph object
#' @param standardized standardize data or not
#' @param nodelabels vector of node names
#' @param alpha significance level for covariate testing
#'
#' @return adjacency matrix of (P)DAG
#' @examples
#' prune_dag(data, curr_amat, alpha)
#' @export
#'
prune_dag <- function(data, curr_amat, alpha=0.00001, standardized=FALSE, nodelabels=colnames(data)){

  start_time = Sys.time()

  # standardize data if needed
  if(standardized==FALSE){
    data <- apply(data, 2, function(x) if(is.numeric(x)){scale(x, center=TRUE, scale=TRUE)} else x)
  }
  data <- as.data.frame(data)

  # create bnlearn object
  curr_graph = empty.graph(colnames(curr_amat))
  amat(curr_graph) = curr_amat
  cond_set <- as.data.frame(matrix(ncol = 0, nrow = 0))
  full_dsep_vec <- c()


  # number of base functions for smoothing function
  p <- dim(as.matrix(data))
  if(p[1]/p[2] < 3*10){
    k_basis <- ceiling(p[1]/(3*p[2]))
  }else{
    k_basis <- min(7, ceiling(p[2]/5))
  }

  gam.control.param = list(nthreads=4, epsilon=1e-04, maxit=100, mgcv.tol=1e-04)

  # try to extend the CPDAG to a DAG and estimate a topological order
  # if no DAG extension is possible, then use random order
  node_order <- nodelabels[sample(c(1:length(nodelabels)))]
  attempt <- 1
  found_dag <- F
  while(found_dag==F && attempt < 6){
    attempt <- attempt + 1
    temp_dag <- tryCatch(bnlearn::pdag2dag(curr_graph, ordering=nodelabels[sample(c(1:length(nodelabels)))]),
                         error=function(e) e, warning=function(w) w)
    if(!is(temp_dag, 'error')){
      found_dag <- TRUE
      node_order <- names(topological.sort(as.igraph(temp_dag)))
    }
  }


  # begin pruning the graph
  curr_struct <- curr_graph$nodes
  for(curr_node in node_order){
    # this include parents and neighbors (undirected edges)
    curr_node_pa <- get_curr_pa(curr_struct, curr_node)
    if(is.null(curr_node_pa)){
      next # no parents/no node
    }
    is_confirmed_pa <- curr_node_pa[,2] # col of true/false
    curr_node_pa <- curr_node_pa[,1] # col of parent/neighbor nodes
    # run gam model
    fm <- get_gam_formula(curr_node, curr_node_pa, k_basis)
    m1 <- gam(fm, method = "REML", optimizer = c("outer","newton"), select = TRUE, data=as.data.frame(data), gam.control=gam.control.param)
    # check if any variable is insignificant, remove edge if so
    m1_pval <- summary(m1)$s.table[,4]
    insig_edges <- which(m1_pval >= alpha)
    insig_edges_pval <- unname(m1_pval[insig_edges])


    # iteratively delete insignificant edges
    # delete from bnlearn object, then update "curr_struct" with new structure
    for(j in 1:length(insig_edges)){
      if(length(insig_edges) == 0){
        break
      }
      insig_node <- curr_node_pa[insig_edges[j]]
      # update skeleton by dropping arc & updating conditioning set, with info on if edge was directed or undirected
      dsep_set <- curr_node_pa[-insig_edges]
      cond_set <- rbind(cond_set, c(insig_node, curr_node, is_confirmed_pa[insig_edges[j]]))
      full_dsep_vec <- c(full_dsep_vec, list(dsep_set))
      # update the bnlearn obj, depending on relation (check if insig node is parent of neighbor)
      if(is_confirmed_pa[insig_edges[j]]){
        curr_graph <- drop.arc(curr_graph, from=insig_node, to=curr_node)
      }
    }
    curr_struct <- curr_graph$nodes
  }
  # rename conditioning set
  if(nrow(cond_set)>0){
    colnames(cond_set) <- c("from", "to", "directed_edge")
    cond_set$dsep <- full_dsep_vec
    # process insignificant undirected edges
    if(sum(cond_set[,3]=="FALSE") > 0){
      cond_set_nbr <- cond_set[cond_set[,3]=="FALSE", ]
      for(i in c(1:nrow(cond_set_nbr))){
        temp_edge <- c(cond_set_nbr[i, 2], cond_set_nbr[i, 1], cond_set_nbr[i, 3], list(cond_set_nbr[i, 4]))
        curr_graph <- drop.arc(curr_graph, from=cond_set_nbr[i, 1], to=cond_set_nbr[i, 2])
        # add opposite direction to cond set
        cond_set <- rbind(cond_set, temp_edge)
      }
    }
    # add everything to blacklist
    bl <- cond_set[,c(1,2)]
    bl <- matrix(bl, ncol = 2, byrow=TRUE, dimnames = list(NULL, c("from", "to")))
    if(is.null(curr_graph$learning$blacklist)){
      curr_graph$learning$blacklist <- bl
    }else{
      curr_graph$learning$blacklist <- rbind(curr_graph$learning$blacklist, bl)
    }
  }

  # record time for GAM procedure
  end_time = Sys.time()
  total_time <- as.numeric(end_time - start_time, unit = "secs")

  return(amat=amat(curr_graph))
}






