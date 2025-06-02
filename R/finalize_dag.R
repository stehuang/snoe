
#' (FOR PRACTICAL PURPOSES) STEP 4 OF ALGORITHM: convert PDAG to DAG
#'
#' In practice, the graph after pruning may contain undirected edges
#' Convert PDAG to DAG by applying likelihood test to determine true causal direction
#'
#' @param data observed data
#' @param curr_amat adjacency matrix of graph object
#' @param alpha significance level for likelihood test
#' @param ll_test_approach "2fold" for running the 2 fold test, other for 5-5 train test split approach
#' @param nodelabels vector containing variable names
#' @param train_ratio percentage of data used for training models
#'
#' @return adjacency matrix of DAG

#' @examples
#' finalize_dag(data, curr_amat=adjmat, ll_test_approach='55')
#' @export
#'
finalize_dag = function(data, curr_amat, nodelabels = colnames(data), alpha = 0.05, train_ratio=0.5, ll_test_approach=NULL){
  # create bnlearn object
  adjmat = curr_amat
  curr_graph = empty.graph(colnames(curr_amat))
  amat(curr_graph) = curr_amat

  # check if this is already a DAG; screen for undirected edges
  if(sum(adjmat + t(adjmat) == 2) == 0){
    out = list("amat"= adjmat, "sigpair" = NULL,
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

  # adding column for storing p-values,diff in ll, and whether edge is oriented
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

  # order edges by number of nbrs + mutual information
  for(curr_n_nbr in c(0:max_common_nbr)){
    udr_edges_subset <- udr_edges[udr_edges[,"n_nbr"]==curr_n_nbr,,drop=F]
    if(nrow(udr_edges_subset)==0) next
    # if there are common nbrs, rank edges
    # if(curr_n_nbr > 0 && nrow(udr_edges_subset > 1)){
    if(nrow(udr_edges_subset > 1)){
      min_mi_vec <- apply(udr_edges_subset, MARGIN=1,
                          FUN=function(x) has_common_pa(dat, adjmat=adjmat, curr_edge=x,
                                                        k=k_basis, compute_min_indep_only = TRUE))
      udr_edges_subset <- cbind(udr_edges_subset, min_mi_vec)
      colnames(udr_edges_subset)[ncol(udr_edges_subset)] <- "min_mi"
      udr_edges_subset <- udr_edges_subset[order(udr_edges_subset[, "min_mi"], decreasing=FALSE),,drop=FALSE]
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
        sig_ind = rbind(sig_ind, c(pa_ind, ch_ind))
        udr_edges[orig_df_index, "has_orient"] <- 1

        # update whitelist in cpdag
        gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, ch_ind)
        # apply Meek's rule (ie convert to cpdag) to orient more edges
        # gam_cpdag <- cpdag(gam_cpdag, wlbl=TRUE)
        adjmat = amat(gam_cpdag)

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

        # update number of neighbors
        udr_edges[,"n_nbr"] <- count_common_nbr(adjmat, udr_edges[,c(1,2)], existing_nbr=udr_edges[,"n_nbr"])
      }else{next}
      gc()
    }
  }


  # second iteration: undirected edges + edges that are sig but couldn't orient
  # run likelihood test again
  udr_edges_subset <- NULL
  if(sum(udr_edges[,"has_eval"] == 0) + sum(udr_edges[,"has_orient"]==0) > 0){
    # edges not evaluated + edge found sig but weren't oriented
    udr_edges_subset <- udr_edges[(udr_edges[,"has_eval"]==0) | (udr_edges[,"has_orient"]==0),,drop=FALSE]
    udr_edges_subset <- udr_edges_subset[order(udr_edges_subset[,"n_nbr"], decreasing=FALSE),,drop=FALSE]
    if(nrow(udr_edges_subset)==0) break
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
        cat("orienting curr root node and edge: ", nodelabels[pa_ind], nodelabels[ch_ind])

        # orient all edges of root node outwards
        pa_node_nbr <- nodelabels[adjmat[,nodelabels[pa_ind]] == adjmat[nodelabels[pa_ind],] & adjmat[,nodelabels[pa_ind]] == 1]
        print(pa_node_nbr)
        if(length(pa_node_nbr) > 0){
          for(nbr_ind in pa_node_nbr){
            if(creates_cycles(adjmat, nbr_ind, pa_ind)==TRUE) next
            gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[pa_ind], to=nbr_ind)
            adjmat = amat(gam_cpdag)
            # update whitelist in cpdag
            gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, pa_ind, which(nodelabels == nbr_ind))
            cat("orienting root node neighbor: ", nodelabels[pa_ind], nbr_ind)
          }
        }

        # # update existing CPDAG and orient edges again
        # gam_cpdag <- set.arc(gam_cpdag, from=nodelabels[pa_ind], to=nodelabels[ch_ind])
        # sig_ind = rbind(sig_ind, c(pa_ind, ch_ind))
        # udr_edges[orig_df_index, "has_orient"] <- 1

        # # update whitelist in cpdag
        # gam_cpdag <- update_bnlearn_dag_wl(gam_cpdag, adjmat=NULL, pa_ind, ch_ind)
        # apply Meek's rule (ie convert to cpdag) to orient more edges
        # gam_cpdag <- cpdag(gam_cpdag, wlbl=TRUE)
        adjmat = amat(gam_cpdag)
        # # update number of neighbors
        # udr_edges[,"n_nbr"] <- count_common_nbr(adjmat, udr_edges[,c(1,2)], existing_nbr=udr_edges[,"n_nbr"])
      }else{next}
      gc()
    }
  }

  # third iteration: only check for undirected edges
  # orient or remove as a result
  if(sum(udr_edges[,"has_orient"]==0)>0){
    # get remaining arcs
    udr_edges_subset <- udr_edges[udr_edges[,"has_orient"]==0,,drop=FALSE]
    udr_edges_subset <- udr_edges_subset[order(udr_edges_subset[,"pval"], decreasing=FALSE),,drop=FALSE]
    if(nrow(udr_edges_subset)==0) break
    for(i in c(1:nrow(udr_edges_subset))){
      pa_ind = udr_edges_subset[i, "pa"]
      ch_ind = udr_edges_subset[i, "ch"]
      orig_df_index <- which(udr_edges[,"pa"] == pa_ind & udr_edges[,"ch"] == ch_ind, arr.ind=TRUE)
      likely_dir <- creates_cycles(adjmat, ch_ind, pa_ind, cpdag=F)
      other_dir <- creates_cycles(adjmat, pa_ind, ch_ind, cpdag=F)
      if(likely_dir == FALSE){
        adjmat[pa_ind, ch_ind] <- 1
        adjmat[ch_ind, pa_ind] <- 0
      }else if(other_dir == FALSE){
        adjmat[pa_ind, ch_ind] <- 0
        adjmat[ch_ind, pa_ind] <- 1
      }else{ # remove undirected edge from graph
        adjmat[pa_ind, ch_ind] <- 0
        adjmat[ch_ind, pa_ind] <- 0
      }
      udr_edges[orig_df_index, "has_orient"] <- 1
    }
  }

  amat(gam_cpdag) <- adjmat

  ll_test_end <- Sys.time()
  if(!is.null(sig_ind)){
    colnames(sig_ind) <- c("from", "to")
    pairs <- nrow(sig_ind)
    sig_ind[,1] <- nodelabels[as.numeric(sig_ind[,1])]
    sig_ind[,2] <- nodelabels[as.numeric(sig_ind[,2])]
  }

  ll_end <- Sys.time()
  total_time <- as.numeric(ll_end-total_start, unit = "secs")
  edgeorient_time <- as.numeric(ll_test_end-ll_test_start, unit = "secs")

  return(adjmat)
}

