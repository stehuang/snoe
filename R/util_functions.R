

# loads rda file
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#' Calculate jacard index (JI) of learned DAG
#'
#' JI: percentage of correct edges among the union of edges in the estimated & true graphs
#'
#' @param true_amat adjacency matrix of true DAG
#' @param temp_amat adjacency matrix of learned DAG
#'
#' @return return the jacard index
#'
#' @examples
#' get_ji(true_dag, learned_dag)
#'
#' @export
get_ji <- function(true_amat, temp_amat){
  true_cpdag = empty.graph(colnames(true_amat))
  amat(true_cpdag) = true_amat

  temp_cpdag = empty.graph(colnames(temp_amat))
  amat(temp_cpdag) = temp_amat

  # get true pos and remove duplicated edges (undirected edges are counted twice)
  tp <- bnlearn::compare(true_cpdag, temp_cpdag, arcs=TRUE)$tp
  tp <- nrow(matrix(tp[!duplicated(t(apply(tp, 1, sort))),], ncol=2))


  # get edge set size of true restricted cpdag & learned cpdag
  true_total_edge <- nrow(directed.arcs(true_cpdag)) + nrow(undirected.arcs(true_cpdag))/2
  temp_total_edge <- nrow(directed.arcs(temp_cpdag)) + nrow(undirected.arcs(temp_cpdag))/2
  ji <- tp/(true_total_edge+temp_total_edge-tp)

  return(ji)
}


#' Calculate results of learned DAG
#'
#' SHD, F1 score, true positives, false positives, false negatives, incorrectly oriented edges
#'
#' @param true_amat ground truth graphical structure adjmat
#' @param temp_amat_list LIST of learned graphical structures adjmat
#'
#' @return vector/matrix containing results; each row corresponds to a learned dag
#'
#' @examples
#' get_f1(true_dag, learned_dag_list)
#'
#' @export
get_f1 <- function(true_amat, temp_amat_list){
  res <- c()
  # convert true amat to bnlearn obj
  true_cpdag = empty.graph(colnames(true_amat))
  amat(true_cpdag) = true_amat

  for(i in c(1:length(temp_amat_list))){
    # convert learned amat to bnlearn obj
    temp_amat <- temp_amat_list[[i]]
    temp_cpdag = empty.graph(colnames(temp_amat))
    amat(temp_cpdag) = temp_amat

    # get true pos and remove duplicated edges (undirected edges are counted twice; should be one edge only)
    tp <- bnlearn::compare(true_cpdag, temp_cpdag, arcs=TRUE)$tp
    tp <- sum(!duplicated(t(apply(tp, 1, sort))))


    # get edge set size of true restricted cpdag & learned cpdag
    true_total_edge <- nrow(directed.arcs(true_cpdag)) + nrow(undirected.arcs(true_cpdag))/2
    # number of undirected edges in graph
    n_undir <- nrow(undirected.arcs(temp_cpdag))/2
    temp_total_edge <- nrow(directed.arcs(temp_cpdag)) + n_undir

    precision <- tp/temp_total_edge
    recall <- tp/true_total_edge
    f1 <- 2*precision*recall/(precision+recall)

    # get incorrectly oriented edges; check w/ skeleton
    tp_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag),
                                    bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$tp
    tp_skeleton <- matrix(tp_skeleton[!duplicated(t(apply(tp_skeleton, 1, sort))),], ncol=2)

    # number of incorrectly oriented edges
    r <- 0
    # number of incorrectly oriented edges that are undirected (supposedly should be directed)
    r_undir <- 0
    if(nrow(tp_skeleton) > 0){
      for(i in 1:nrow(tp_skeleton)){
        n1 <- tp_skeleton[i,1]
        n2 <- tp_skeleton[i,2]
        # need edge direction to be same
        if(true_amat[n1,n2]!=temp_amat[n1,n2] || true_amat[n2,n1]!=temp_amat[n2,n1]){
          r <- r+1
          # check if wrong dir is b/c undirected edge
          if(temp_amat[n1,n2] + temp_amat[n2,n1]==2){
            r_undir <- r_undir + 1
          }
        }
      }
    }
    fp <- temp_total_edge-tp-r
    fn <- true_total_edge-tp-r

    # another way to calculate f1
    f1_2 <- 2*tp/(2*tp+fp+fn+2*r)

    # calculate SHD manually; fp+fn+wrong_dir
    shd_score <- fp+fn+r

    res <- rbind(res, c(shd_score, f1, tp, fp, fn, r, r_undir, n_undir))
    colnames(res) <- c("shd","f1", "tp", "fp", "fn", "r", "r_undir", "n_undir")
  }
  return(res)
}


# calculate number of incorrectly oriented edges
get_wrong_dir <- function(true_cpdag, temp_cpdag){
  true_amat <- amat(true_cpdag)
  temp_amat <- amat(temp_cpdag)
  # get incorrectly oriented edges
  tp_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag),
                                  bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$tp
  tp_skeleton <- matrix(tp_skeleton[!duplicated(t(apply(tp_skeleton, 1, sort))),], ncol=2)

  r_edges <- c()
  r <- 0
  num_temp_undir <- 0
  num_temp_undir2 <- 0
  for(i in 1:nrow(tp_skeleton)){
    if(nrow(tp_skeleton)==0){
      break
    }
    n1 <- tp_skeleton[i,1]
    n2 <- tp_skeleton[i,2]
    same_orientation <- (true_amat[n1,n2]==temp_amat[n1,n2]) + (true_amat[n2,n1]==temp_amat[n2,n1])
    if(same_orientation != 2){
      r <- r+1
      if(temp_amat[n1,n2]==temp_amat[n2,n1]){ # check if edge in learned DAG is undirected
        num_temp_undir2 <- num_temp_undir2+1
      }
    }
    if(true_amat[n1,n2]!=true_amat[n2,n1]){
      num_temp_undir <- num_temp_undir+1
    }
  }

  # get wrong dir, and directed edges
  return(c(r, r-num_temp_undir2, num_temp_undir2))
}


# returns: percentage of undirected edges captured
get_undir_pct <- function(true_amat, cpdag_amat, temp_amat){
  # number of undirected edges undetected
  N.det.edge = sum((abs(temp_amat - cpdag_amat) + abs(true_amat - cpdag_amat))==2)
  # number of undirected edges total
  N.undir.edge = nrow(matrix(which((cpdag_amat == t(cpdag_amat) & (cpdag_amat != 0)), arr.ind = TRUE), ncol = 2))/2
  # percentage of undirected edges from CPDAG captured by GAM
  return(N.det.edge/N.undir.edge)
}


# computes shd, f1, ji, and number of nonlinear edgs missing
analyze_cpdag <- function(temp.cpdag, true.cpdag){
  shd_cpdag <- bnlearn::shd(true=true.cpdag, learned=temp.cpdag, wlbl=TRUE)
  f1_cpdag <- get_f1(true.cpdag, temp.cpdag)
  ji_cpdag <- get_ji(true.cpdag, temp.cpdag)
  noninv_fn <- count_fn(true.cpdag, temp.cpdag)
  # return(c(shd_cpdag, f1_cpdag, ji_cpdag, noninv_fn, temp$gam_edge_time, temp$edgeorient_time))
  return(c(shd_cpdag, f1_cpdag, ji_cpdag, noninv_fn))
}


# get number of missing edges that are non-invertible/non-linear
count_fn <- function(true_cpdag, temp_cpdag){
  true_amat <- amat(true_cpdag)
  temp_amat <- amat(temp_cpdag)

  fn_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag),
                                  bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$fn
  if(!is.null(true_cpdag$learning$whitelist)){
    non_inv_fn <- nrow(inner_join(as.data.frame(true_cpdag$learning$whitelist), as.data.frame(fn_skeleton), by = c("from", "to")))
  }else{
    non_inv_fn <- 0
  }
  return(non_inv_fn)
}


#' Imports pre-saved/bnlearn graphs as bnlearn object
#'
#' Includes self-defined graphs and graphs contained in .rda files
#'
#' @param network name of network
#' @param fixed_edges true/statement on whether to keep certain edges fixed in cpdag
#'
#' @return vector/matrix containing results; each row corresponds to a learned dag
#'
#' @examples
#' get_true_dag('mildew')
#'
#' @export
get_true_dag <- function(network='child', fixed_edges=TRUE){
  dag_obj <- NULL
  avail_bn <- c(
    # Discrete Bayesian networks
    "asia", "cancer", "earthquake", "sachs", "survey",
    "alarm", "barley", "child", "insurance", "mildew", "water",
    "hailfinder", "hepar2", "win95pts",
    "munin", "munin1", "munin2", "munin3", "munin4",
    "andes", "diabetes", "link", "pathfinder", "pigs",
    # Gaussian Bayesian networks
    "ecoli70", "magic-niab", "magic-irri", "arth150",
    # Conditional linear Gaussian Bayesian networks
    "healthcare", "sangiovese", "mehra"
  )
  if (!network %in% avail_bn){

    stop(sprintf("network must be one of %s",
                 paste(avail_bn, collapse = ", ")))
  }

  x1 <- ifelse(network %in% sprintf("munin%s", seq_len(4)), "munin4", network)
  x2 <- ifelse(network == "mehra", "mehra-complete", network)

  dag_obj <- readRDS(
    file = url(sprintf("https://www.bnlearn.com/bnrepository/%s/%s.rds",
                       x1, x2))
  )
  true_amat <- amat(dag_obj)
  true.dag = empty.graph(colnames(true_amat))
  amat(true.dag) = true_amat
  labels = paste("X", 1:length(true.dag$nodes), sep = "")
  true.dag <- rename.nodes(true.dag, labels)

  return(true.dag)
  # if(network=="1"){
  #   # EXAMPLE 1
  #   true_amat <- rbind(c(0,1,1,0), c(0,0,0,0), c(0,0,0,1), c(0,0,0,0))
  #   colnames(true_amat) <- c("x1","x2","x3","x4")
  #   rownames(true_amat) <- c("x1","x2","x3","x4")
  #   ninv.arcs <- rbind(c("x1", "x2"), c("x1", "x3"))
  # }else if(network=="2"){
  #   # EXAMPLE 2
  #   # same network as network 3, but with node 6 removed
  #   true_amat <- rbind(c(0,1,0,0,0), c(0,0,1,1,0),c(0,0,0,0,0),
  #                      c(0,0,0,0,0), c(0,0,0,1,0))
  #   colnames(true_amat) <- c("x1","x2","x3","x4","x5")
  #   rownames(true_amat) <- c("x1","x2","x3","x4","x5")
  #   # ninv.arcs <- rbind(c("x2", "x3"))
  #   ninv.arcs <- c("x2", "x4")
  # }else if(network=="3"){
  #   # EXAMPLE 3
  #   # true dag: x1->x2->x3, x2->x4<-x5, x4->x6
  #   true_amat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0),
  #                      c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
  #   colnames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
  #   rownames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
  #   ninv.arcs <- rbind(c("x2", "x4"))
  # }else if(network=="asia"){
  #   dag_obj <- loadRData("asia.rda")
  #   ninv.arcs <- rbind(c("tub", "either"),c("smoke","lung"))
  # }else if(network=="sachs"){
  #   dag_obj <- loadRData("/Users/StellaHuang/Documents/UCLA/Research/GAM/sachs.rda")
  #   ninv.arcs <- rbind(c("PKC", "Mek"),c("PIP3","PIP2"), c("PKA","P38"), c("PKA","Erk"))
  # }else if(network=="child"){
  #   dag_obj <- loadRData("child.rda")
  #   ninv.arcs <- rbind(c("BirthAsphyxia", "Disease"),c("LungParench","ChestXray"),
  #                      c("CardiacMixing","HypoxiaInO2"), c("Disease","Sick"))
  # }else if(network=="sang"){
  #   dag_obj <- loadRData("sangiovese.rda")
  # }else if(network=="mehra"){
  #   dag_obj <- loadRData("mehra.rda")
  # }else if(network=="mildew"){
  #   dag_obj <- loadRData("mildew.rda")
  #   true_amat <- amat(mildew)
  #   # c("Disease","LungParench") c("BirthAsphyxia", "Disease")
  # }else if(network=="alarm"){
  #   dag_obj <- loadRData("alarm.rda")
  #   true_amat <- amat(insurance)
  # }else if(network=="water"){
  #   dag_obj <- loadRData("water.rda")
  #   true_amat <- amat(water)
  # }else if(network=="magic2"){
  #   dag_obj <- loadRData("magic2.rda")
  # }else if(network=="ecoli"){
  #   dag_obj <- loadRData("ecoli.rda")
  # }else if(network=="hf"){
  #   dag_obj <- loadRData("hailfinder.rda")
  # }else if(network=="magic"){
  #   dag_obj <- loadRData("magic.rda")
  # }else if(network=="hepar2"){
  #   dag_obj <- loadRData("hepar2.rda")
  # }else if(network=="pathfinder"){
  #   dag_obj <- loadRData("pathfinder.rda")
  # }else if(network=="win"){
  #   dag_obj <- loadRData("win.rda")
  # }else if(network=="water_rand_order"){
  #   dag_obj <- loadRData("water.rda")
  # }else if(network=="magic2_rand_order"){
  #   dag_obj <- loadRData("magic2.rda")
  # }else if(network=="mehra_rand_order"){
  #   dag_obj <- loadRData("mehra.rda")
  # }else if(network=="magic_rand_order"){
  #   dag_obj <- loadRData("magic.rda")
  # }else if(network=="mildew_rand_order"){
  #   dag_obj <- loadRData("mildew.rda")
  # }else if(network=="alarm_rand_order"){
  #   dag_obj <- loadRData("alarm.rda")
  # }else{
  #   print("else network 3 chosen; creating true amat")
  #   true_amat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0),
  #                      c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
  #   colnames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
  #   rownames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
  # }
  #

  # create bnlearn obj of true dag
  # if(!is.null(dag_obj)) true_amat <- amat(dag_obj)
  # true.dag = empty.graph(colnames(true_amat))
  # amat(true.dag) = true_amat
  # nodeN = length(true.dag$nodes)
  # labels = paste("X", 1:nodeN, sep = "")
  # true.dag <- rename.nodes(true.dag, labels)
  #
  # return(true.dag)
}


# obsolete function
# for classifying edges in learned dag
classify_edges <- function(init_cpdag, final_cpdag, nonlinear_edges){
  init_edges <- cbind(init_cpdag$arcs, "nonlinear"= rep(0,nrow(init_cpdag$arcs)))
  final_edges <- final_cpdag$arcs
  # get linear edges (i.e. from initial learning alg and preserved till the end)
  linear_edges <- bnlearn::compare(final_cpdag, init_cpdag,  arcs=TRUE)$tp
  linear_edges <- cbind(linear_edges, "nonlinear"= rep(0,nrow(linear_edges)))
  # get nonlinear edges (i.e. found via oriented, added procedures)
  if(!is.null(nonlinear_edges)){
    nonlinear_edges <- cbind(nonlinear_edges, "nonlinear"= rep(1,nrow(nonlinear_edges)))

    # left join to get linear/nonlinear status
    final_edges <- merge(x=final_edges,y=rbind(linear_edges, nonlinear_edges),
                         by=c("from","to"), all.x=TRUE)
    final_edges <- final_edges[!duplicated(final_edges[,c(1,2)], fromLast=TRUE),]
    # if edge status is missing, then it's linear
    final_edges[is.na(final_edges)] <- 0
  }else{
    final_edges <- cbind(final_edges, "nonlinear"= rep(0,nrow(final_edges)))
  }

  return(final_edges)
}



#' Extends CPDAG to DAG
#'

#' @param cpdag bnlearn object of the cpdag
#'
#' @return bnlearn object that is a DAG; if no extension if found, then original cpdag is returned
#'
#' @examples
#' cpdag2dag(curr_cpdag)
#'
#' @export
cpdag2dag <- function(cpdag){
  # convert CPDAG to random DAG in equivalence class
  counter <- 0
  # need to supply ordering, so we pass in random order
  # permutate order until cpdag is converted into dag
  while(counter < 100){
    dag <- try(bnlearn::pdag2dag(cpdag, names(cpdag$nodes)[sample(1:length(names(cpdag$nodes)))]), silent=TRUE)
    if(class(dag) != "try-error"){
      return(dag)
    }
    counter <- counter+1
  }
  return(cpdag)
}




