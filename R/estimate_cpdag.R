


# function: get local structure that corresponds to the "nodes" feature in bnlearn obj
# input: adjacency matrix
# output: local structure + dsep_set
get_local_structure <- function(adjmat, add_dsep_set = TRUE){
  nodelabels <- colnames(adjmat)
  # get arcs in skeleton
  skel_arcs <- matrix(which(adjmat+t(adjmat)>0, arr.ind = T), ncol=2)
  skel_arcs[,1] <- rownames(adjmat)[as.numeric(skel_arcs[,1])]
  skel_arcs[,2] <- colnames(adjmat)[as.numeric(skel_arcs[,2])]
  skel_arcs <- skel_arcs[skel_arcs[,1]!=skel_arcs[,2],]
  colnames(skel_arcs) <- c("from", "to")
  skeleton <- bnlearn:::cache.structure(nodelabels, skel_arcs)
  # save separating set
  if(add_dsep_set){
    dsep_set <- get_dsep_set(adjmat)
    attr(skeleton, "dsep.set") <- dsep_set
  }
  return(skeleton)
}


# function: add "dsep.set" attribute to estimated skeleton; used for finding new separating sets
# input: adjmat of skeleton
# output: list of lists to be added as "dsep.set" attribute
get_dsep_set <- function(adjmat, nodelabels=rownames(adjmat)){
  # store dsep_set results
  dsep.set <- list()
  # get deleted, undirected edges
  # skel_arcs <- matrix(which(adjmat==0, arr.ind = T), ncol=2)
  skel_arcs <- matrix(which(adjmat+t(adjmat)==0, arr.ind = T), ncol=2)
  skel_arcs[,1] <- rownames(adjmat)[as.numeric(skel_arcs[,1])]
  skel_arcs[,2] <- colnames(adjmat)[as.numeric(skel_arcs[,2])]
  skel_arcs <- skel_arcs[!duplicated(t(apply(skel_arcs, 1, sort))), ] # remove duplicates from undirected edges
  skel_arcs <- skel_arcs[skel_arcs[,1]!=skel_arcs[,2],]
  colnames(skel_arcs) <- c("from", "to")

  # find the smallest separating set (i.e. smallest number of neighboring nodes) for two nodes connected by edge
  for(i in c(1:nrow(skel_arcs))){
    arc <- as.vector(skel_arcs[i,])
    n1 <- arc[1]
    n2 <- arc[2]
    sep_set <- as.vector(nodelabels[adjmat[,n1]==1])
    sep_set2 <- as.vector(nodelabels[adjmat[,n2]==1])
    if(adjmat[n1,n2] + adjmat[n2,n1]==2){ # undirected edge; consider both sep set
      if(length(sep_set) > length(sep_set2)){ # take the smallest separating set, including an empty one
        sep_set <- sep_set2
      }
    }else{
      if(adjmat[n1,n2]==0){ # if n1 is not a parent of n2, then use sep set of n2
        sep_set <- sep_set2
      }
    }
    # p-value is needed for formatting reasons in vstruct detection; put placeholder here
    dsep.set <- append(dsep.set, list(list(arc=c(n2, n1), p.value=0.5, dsep.set=sep_set)))
  }
  return(dsep.set)
}


#' STEP 1 OF ALGORITHM: Learning the initial CPDAG from data
#'
#' Learns a dense CPDAG from data; skeleton is estimated under a relaxed threshold, and conditional indep relations are detected using a strigent threshold
#'
#' @param data observed data
#' @param learning_alg options are "PC" or "ccdr"
#' @param pc_params parameters to pass into PC algorithm; alpha is threshold for CI tests for v-structure detection, alpha_dense is for learning the skeleton
#' @param ccdr_pc_alpha estimates CPDAG from PC to determine size of DAG from ccdr
#' @param nodelabels vector of variable names
#' @param cluster use cluster for computing or not (bnlearn function argument)
#' @param debug display debug process or not (bnlearn function argument)
#' @param ccdr_params parameters to pass into ccdr algorithm
#'
#' @return adjacency matrix of learned cpdag
#' @examples
#' estimate_cpdag(data, standardized=T, learning_alg="PC", pc_params(alpha=0.05, alpha_dense=0.25, test='cor'))
#' @export
#'
#'

estimate_cpdag <- function(data, nodelabels=colnames(data), cluster=NULL,
                           learning_alg='PC', pc_params=list(alpha=0.01, alpha_dense=0.25, test='cor'), ccdr_pc_alpha = 0.01,
                           ccdr_params=list(pc_alpha=ccdr_pc_alpha, lambdas.length=15, gamma=2, max.iters = sparsebnUtils::default_max_iters(ncol(data)),
                                            ccdr_alpha = sparsebnUtils::default_alpha(), error.tol = 1e-4, verbose=FALSE), debug=FALSE){
  start_time = Sys.time()

  # standardize data if needed
  data <- as.data.frame(data)
  data <- bnlearn:::check.data(data, allow.missing = TRUE)

  # run PC algorithm to estimate skeleton
  if(learning_alg %in% c("PC", "pc")){ ## pc algorithm
    pc_alpha <- pc_params$alpha
    pc_alpha_dense <- pc_params$alpha_dense
    test <- pc_params$test
    start_time_est = Sys.time()
    # discretize data if using mutual information
    if(test=="mi"){
      data <- apply(data, MARGIN=2, FUN=function(x) as.factor(infotheo::discretize(x, disc="equalwidth", nbins=round(nrow(data))/50)$X))
      data <- as.data.frame(unclass(data),stringsAsFactors=TRUE)
      data <- bnlearn:::check.data(data, allow.missing = TRUE)
    }
    dense_skeleton <- bnlearn:::pc.stable.backend(x=as.data.frame(data), whitelist = NULL, blacklist = NULL, test=test,
                                                  alpha=pc_alpha_dense, max.sx = 3, B=0L)
    # remove dsep set to re-run CI test with stricter threshold
    attr(dense_skeleton, "dsep.set") <- NULL
    # detect new v-structures with strigent threshold
    dense_skeleton_obj <- learn_arc_directions_pdag(x=as.data.frame(data), local.structure = dense_skeleton,
                                                    blacklist=NULL, whitelist=NULL, max.sx = 3, debug=debug,
                                                    test=test, alpha=pc_alpha, vs=NULL)
    pdag_amat <- bnlearn:::arcs2amat(dense_skeleton_obj$arcs, colnames(data))
    init_cpdag <- empty.graph(nodelabels)
    amat(init_cpdag, check.cycles=FALSE) <- pdag_amat

    # add directed edges to whitelist (ensures edge is strongly protected)
    init_cpdag$learning$whitelist <- bnlearn::directed.arcs(init_cpdag)
    end_time_est = Sys.time()
  }else{ # run sparsebn to obtain initial cpdag
    start_time_est = Sys.time()
    # estimate solution path; sb_dag returns multiple DAGs
    dat <- sparsebnData(data, type = "c", ivn = NULL)
    sb_dag <- ccdrAlgorithm::ccdr.run(data = dat,
                                      lambdas.length = ccdr_params$lambdas.length,
                                      whitelist = NULL,
                                      blacklist = NULL,
                                      gamma = ccdr_params$gamma,
                                      error.tol = ccdr_params$error.tol,
                                      max.iters = ccdr_params$max.iters,
                                      alpha = ccdr_params$ccdr_alpha,
                                      verbose = ccdr_params$verbose)
    dag_index <- sparsebnUtils::select.parameter(sb_dag, dat)
    # select optimal dag based on PC skeleton density
    pcnet <- pc.stable(data, alpha=ccdr_pc_alpha)
    pcnet_edges <- nrow(bnlearn::skeleton(pcnet)$arcs)/2
    # select dag with the closest number of edges as pc skeleton
    curr_dag <- sparsebnUtils::select(sb_dag, edges = pcnet_edges*(1))
    init_skeleton <- sparsebnUtils:::to_bn(curr_dag)$edges
    init_cpdag <- get_init_cpdag(data, skeleton=init_skeleton, blacklist=NULL, debug=FALSE)
    end_time_est = Sys.time()
  }
  end_time = Sys.time()
  total_time <- as.numeric(end_time - start_time, unit = "secs")
  return(amat = amat(init_cpdag))
}



# function: extract colliders nodes in skeleton; modified from bnlearn:::colliders.backend
# input: adjacency matrix
# output: local structure (skeleton) + separating set
colliders_backend_mod <- function(x, return.arcs = FALSE, including.shielded = TRUE,
                                  including.unshielded = TRUE, debug = FALSE) {

  nodes = names(x$nodes)

  coll = .Call(bnlearn:::call_colliders,
               arcs = x$arcs,
               nodes = nodes,
               return.arcs = return.arcs,
               shielded = including.shielded,
               unshielded = including.unshielded,
               debug = debug)

  if (return.arcs) {

    coll = bnlearn:::arcs.rbind(coll[, c("X", "Z")], coll[, c("Y", "Z")])
    coll = bnlearn:::arcs.unique(coll, nodes = nodes)

  }#THEN

  return(coll)

}


# function: performs v-structure detection & edge orientation to get PDAG
# input: data, skeleton
# output: adjacency matrix of learned pdag
detect_v_struct <- function(data, skeleton, blacklist=NULL,
                            whitelist=NULL, test='cor', debug=FALSE,alpha=0.05){
  # this function does v-structure detection first, then edge orientation
  cpdag <- bnlearn:::learn.arc.directions(x=data, local.structure = skeleton,
                                          blacklist=blacklist, whitelist=whitelist,
                                          test=test, alpha=alpha, max.sx = NULL, debug=debug)
  # jireh's implementation
  # pdag_ppc <- (x=data, local.structure=skeleton, complete=bnlearn:::check.data(data)$complete.nodes,
  # score='bic-g', blacklist=blacklist, whitelist=whitelist, debug = debug)
  nodes <- names(data)
  # convert arcs to adjacency matrix
  Amat <- bnlearn:::arcs2amat(cpdag$arcs, nodes)
  # Amat <- bnlearn:::arcs2amat(pdag_ppc$arcs, nodes)
  return(Amat)
}


# function: converts DAG to CPDAG with option to exclude edges (and certain orientations)
# input: estimated DAG/skeleton, blacklist of edges
# output: CPDAG of the input structure
get_init_cpdag <- function(data, skeleton, blacklist=NULL, debug=FALSE){
  if(is.null(blacklist) == FALSE){
    blacklist <- as.matrix(blacklist[,c(1,2)])
    colnames(blacklist) <- c("from", "to")
  }

  gam_cpdag = bnlearn::empty.graph(colnames(data))
  # make sure deleted edges aren't added back
  gam_cpdag$learning$blacklist <- blacklist

  if(class(skeleton)!="bn"){ # output from ppC alg, which are edges
    amat(gam_cpdag) = detect_v_struct(data=as.data.frame(data), skeleton=skeleton, blacklist=blacklist,debug=debug)
  }else{
    gam_cpdag <- cpdag(skeleton)
  }
  return(gam_cpdag)
}

# function: converts a PDAG to DAG
# input: pdag, name of nodes
# output: DAG if successful, PDAG if not

#' Converts a PDAG to DAG
#'
#' @param pdag bnlearn object of pdag
#' @param nodelabels vector of variables names
#'
#' @return adjacency matrix of dag if successful, adjacency matrix of pdag otherwise
#' @export
#'
#' @examples
#' extend_to_dag(example_pdag, c("X1", "X2", "X3"))
extend_to_dag <- function(pdag, nodelabels){
  pdag_amat <- bnlearn:::arcs2amat(pdag$arcs, nodelabels)
  temp_cpdag <- empty.graph(nodelabels)
  amat(temp_cpdag, check.cycles=FALSE) <- pdag_amat
  res <- tryCatch(bnlearn::cextend(temp_cpdag), error=function(e) e, warning=function(w) w)
  if(is(res, 'error')){
    print("cannot find extension of PDAG to DAG")
    return(pdag_amat)
  }else{
    return(amat(res))
  }
}



# function: converts skeleton to PDAG with detected v-structures; Meek's rule is not applied yet
# modified from bnlearn:::learn.arc.directions
# input: estimated skeleton, blacklist of edges
# output: PDAG of the input structure
learn_arc_directions_pdag <- function (x, cluster = NULL, local.structure, whitelist, blacklist,
                                       test, alpha, B=NULL, max.sx = ncol(data), debug = FALSE, vs=NULL, return_vs_only=F, n_max=50)
{
  nodes = names(x)
  arcs = bnlearn:::mb2arcs(local.structure, nodes)
  to.drop = !apply(arcs, 1, function(x) {
    bnlearn:::is.blacklisted(blacklist, x)
  })
  arcs = arcs[to.drop, , drop = FALSE]
  if(is.null(vs)){
    vs = do.call("rbind", bnlearn:::vstruct.detect(nodes = nodes, arcs = arcs,
                                                   mb = local.structure, data = x, alpha = alpha,
                                                   test = test, blacklist = blacklist, max.sx = max.sx,
                                                   debug = debug))
  }
  if(return_vs_only) return(vs)
  rownames(vs) = NULL
  if (!is.null(vs)) {
    vs = vs[order(vs[, "max_a"], decreasing = FALSE), ]
    arcs = bnlearn:::vstruct.apply(arcs = arcs, vs = vs, nodes = nodes,
                                   debug = debug)
  }
  learning = list(whitelist = whitelist, blacklist = blacklist,
                  test = test, args = list(alpha = alpha), ntests = test.counter())
  if (!is.null(B)) learning$args$B = B
  pdag = list(learning = learning, nodes = structure(rep(0, length(nodes)), names = nodes), arcs = arcs, vs=vs)
  # orient edges
  pdag_amat <- bnlearn:::arcs2amat(pdag$arcs, colnames(data))
  pdag_amat <- phsl:::apply_cpdag_rules(pdag=pdag_amat, nodes=nodes, remove_invalid=TRUE, debug=debug)


  # check for cycles
  contains_cycles <- bnlearn:::is.acyclic(bnlearn:::amat2arcs(pdag_amat, nodes), nodes, directed=T)
  counter <- 1
  while(contains_cycles==FALSE && counter < 25){
    counter <- counter+1
    vs = vs[sample(c(1:nrow(vs))), ]
    # vs = do.call("rbind", bnlearn:::vstruct.detect(nodes = nodes, arcs = arcs,
    #                                                mb = local.structure, data = x, alpha = alpha,
    #                                                test = test, blacklist = blacklist, max.sx = max.sx,
    #                                                debug = debug))
    arcs = bnlearn:::vstruct.apply(arcs = arcs, vs = vs, nodes = nodes,
                                   debug = debug)
    pdag = list(learning = learning, nodes = structure(rep(0, length(nodes)), names = nodes), arcs = arcs, vs=vs)
    # orient edges
    # pdag <- bnlearn:::cpdag.backend(pdag, fix=T, debug = debug)
    pdag_amat <- bnlearn:::arcs2amat(pdag$arcs, nodes)
    pdag_amat <- phsl:::apply_cpdag_rules(pdag=pdag_amat, nodes=nodes, remove_invalid=TRUE, debug=debug)
    if (!is.null(B)) learning$args$B = B
    contains_cycles <- bnlearn:::is.acyclic(bnlearn:::amat2arcs(pdag_amat, nodes), nodes, directed=T)
    # can_extend_to_dag <- extend_to_dag(pdag, nodes)
  }

  return(pdag)
}


# function: detects v-structures in skeleton using separating set
# modified from bnlearn:::vstruct.detect
# output: PDAG of the input structure

vstruct_detect_mod <- function (nodes, arcs, mb, data, alpha, B = NULL, test, blacklist,
                                max.sx = ncol(data), debug = FALSE)
{
  vstruct.centered.on = function(x, mb, data, dsep.set) {
    if (debug) {
      cat("----------------------------------------------------------------\n")
      cat("* v-structures centered on", x, ".\n")
    }
    tos = arcs[(arcs[, "to"] == x), "from"]
    if (length(tos) < 2)
      return(NULL)
    tos.combs = bnlearn:::subsets(tos, 2)
    vs = NULL
    for (j in 1:nrow(tos.combs)) {
      y = tos.combs[j, 1]
      z = tos.combs[j, 2]
      if (debug)
        cat("  * checking", y, "->", x, "<-", z, "\n")
      if (is.listed(arcs, c(y, z), either = TRUE))
        next
      if (is.null(blacklist))
        blacklisted.arc = FALSE
      else blacklisted.arc = is.listed(blacklist, c(y,
                                                    z), both = TRUE)
      if (!is.null(dsep.set) && !blacklisted.arc) {
        el = dsep.set[[which(sapply(dsep.set, function(x) setequal(x$arc,
                                                                   c(y, z))))]]
        if (x %!in% el$dsep.set) {
          if (debug)
            cat("    @ detected v-structure", y, "->",
                x, "<-", z, "from d-separating set.\n")
          vs = rbind(vs, data.frame(max_a = el$p.value,
                                    y, x, z))
        }
      }
      else {
        sx = smaller(setdiff(mb[[y]][["mb"]], c(x, z)),
                     setdiff(mb[[z]][["mb"]], c(x, y)))
        if (debug)
          cat("    > chosen d-separating set: '", sx,
              "'\n")
        a = allsubs.test(x = y, y = z, fixed = x, sx = sx,
                         data = data, test = test, B = B, alpha = alpha,
                         max = min(max.sx, length(sx)), debug = debug)
        if (a["p.value"] <= alpha) {
          if (debug)
            cat("    @ detected v-structure", y, "->",
                x, "<-", z, "\n")
          vs = rbind(vs, data.frame(max_a = a["max.p.value"],
                                    y, x, z))
        }
      }
    }
    return(vs)
  }
  sapply(nodes, vstruct.centered.on, mb = mb, data = data,
         dsep.set = attr(mb, "dsep.set"), simplify = FALSE)
}


# function: adding the detecting v-structures to the skeleton
# modified from bnlearn:::vstruct.apply
# input: detected edges, v-structures
# output: PDAG of the input structure
vstruct_apply_mod <- function (arcs, vs, nodes, debug = FALSE)
{
  if (debug)
    cat("----------------------------------------------------------------\n")
  for (i in seq(nrow(vs))) {
    x = vs[i, "x"]
    y = vs[i, "y"]
    z = vs[i, "z"]
    max_a = vs[i, "max_a"]
    # set edge direction
    temp = bnlearn:::set.arc.direction(y, x, arcs, debug=debug)
    temp = bnlearn:::set.arc.direction(z, x, temp, debug=debug)
    # do not set if it induces a cycle
    if (!bnlearn:::is.acyclic(temp, nodes, directed = TRUE)) {
      if (debug)
        cat("* not applying v-structure", y, "->", x,
            "<-", z, "(", max_a, ")\n")
      warning("vstructure ", y, " -> ", x, " <- ", z,
              " is not applicable, ", "because one or both arcs introduce cycles in the graph.")
      next
    }
    if (debug)
      cat("* applying v-structure", y, "->", x, "<-",
          z, "(", max_a, ")\n")
    arcs = temp
  }
  return(arcs)
}




