


# function: generates data for a pair (x,y) where y=f(x)
# input: sample size, variable x, error term (T/F), mean & sd of x, function types to use
# output: matrix containing values for (x,y)
pairdatagen = function(n = NULL, x = NULL, error = T, mx = NULL, sdx = NULL, data_fx="inv", ninv.method = "gp", inv.method = "random",
                       se = 1, lo = -1, hi = 1, tau.p = NULL, noninv_split=T, error.distr = 'normal'){
  scale01 = function(x){
    (x-min(x))/diff(range(x))
  }

  if(is.null(mx)) mx = 0 # mean of x
  if(is.null(sdx)) sdx = 1 # sd of x
  if(is.null(x)) x = rnorm(n, mx, sdx)

  ## random error
  if(error){
    se_normal <- runif(1, 0.5, 0.75)
    if(is.null(error.distr)){
      e = rnorm(n, 0, se_normal)
    }else if(error.distr == "normal"){
      e = rnorm(n, 0, se_normal)
    }else if(error.distr == "laplace"){
      e = rlaplace(n=n, location = 0, scale = 1)
    }else if(error.distr == "gumbel"){
      e = rgumbel(n=n, location = 0, scale = 1)
    }else if(error.distr == "t"){
      e = rt(n=n, df=5)
    }
  }else{
    e <- rep(0, n)
  }

  sign = sample(c(-1,1), 1)
  y <- NULL

  if(data_fx=='linear'){
    # simple linear regression
    a = runif(1, lo, hi)
    b = max(runif(10, lo, hi))
    y = a + b*x + e
  }

  if(data_fx=="nonlinear"){
    data_fx <- sample(c("inv", "ninv"), 1)
  }

  if(data_fx=="inv"){
    if(inv.method=="random") inv.method = sample(c("piecewise", "exp", "cubic", "sinh"), 1)
    if(inv.method == "piecewise") {
      tau.t = mean(sample(x, 10))
      x1 = as.matrix(x[which(x < tau.t)])
      x2 = as.matrix(x[which(x >= tau.t)])

      b12 = sample(tail(abs(runif(15, lo, hi))), 2)
      b1 = b12[1] * sign
      b2 = b12[2] * sign
      a1 = runif(1, lo, hi)
      a2 = a1 + tau.t * (b1 - b2)

      y = e
      y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
      y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
    }else if(inv.method == "exp"){
      y = runif(1, 0.05, 0.35)*exp(x) + e
    }else if(inv.method == "cubic"){
      y = sign * x^3 * runif(1, 0.1, 0.5) + e
    }else if(inv.method == "sinh"){
      y = sign * sinh(x) + e
    }else if (inv.method == "tanh"){
      y = sign * tanh(x) + e
    }else if(inv.method == "log"){
      x = scale01(x)+0.01
      y = log(x) + e
    }else if(inv.method == "sigmoid"){
      a <- rexp(1, rate=4)+1
      b <- runif(1, 0.5, 2) * sample(c(-1,1), 1)
      c <- runif(1, -2, 2)
      y <- a * (b*(x+c))/(1+abs(b*(x+c))) + e
    }

    if(is.null(y)) stop("Specified function is not found amongst all invertible functions")
  }

  if(data_fx=="ninv"){
    if(ninv.method == "random") ninv.method = sample(c("piecewise", "quadratic", "cubic", "x4", "gp"), 1)

    # gaussian process
    if(ninv.method == "gp"){
      l = runif(1, 5, 5.25)
      d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
      Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
      # generate 1 sample; returns 500*1 dim vector, since x is 500*1
      y = as.vector(mvtnorm::rmvnorm(n=1, sigma=Sigma_SE) + e)
    }

    if(ninv.method == "piecewise") {
      tau.t = unname(quantile(x,probs=c(0.15, 0.2, 0.25, 0.3, 0.7, 0.75, 0.8, 0.85)))
      x1 = as.matrix(x[which(x < tau.t)])
      x2 = as.matrix(x[which(x >= tau.t)])

      b12 = sample(tail(abs(runif(15, lo, hi))), 2)
      b1 = b12[1] * sign
      b2 = b12[2] * sign * -1/4
      a1 = runif(1, lo, hi)
      a2 = a1 + tau.t * (b1 - b2)

      y = e
      y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
      y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
    }


    if(ninv.method == "nl1") y = sample(c(-1,1), 1)*max(runif(2, lo, hi))*x^2 + sample(c(-1,1), 1)*min(runif(2, lo, hi))*x + e

    if(ninv.method == "nl2") y = sample(c(-1,1), 1)*cos(max(2, 2*runif(1, lo, hi))*x) + e

    if(ninv.method == "nl3") y = sample(c(-1,1), 1)*runif(1, lo, hi)*x^3 + 2*runif(1, lo, hi)*x^2 + e

    if(ninv.method == "nl4") y = sample(c(-1,1), 1)*tanh(x) + sample(c(-1,1), 1)*cos(2.5*x) + sample(c(-1,1), 1)*x^2 + e

    if(ninv.method == "cubic") y = sample(c(-1,1), 1)*(x-1)^2*(runif(1, 0.75, 1.25)*x+3) + e

    if(ninv.method == "quadratic") y = sample(c(-1,1), 1)*(x+runif(1, 0.1, 0.65))^2 + e

    if(ninv.method == "x4") y = sample(c(-1,1), 1)*(0.75*x^4-2*x^2-x^3) + e

    if(is.null(y)) stop("Specified function is not found amongst all non-invertible functions")
  }

  dat = NULL
  dat$x = scale(x) #scale(x)
  dat$y = scale(y) #scale(y)
  return(dat)
}


# function: generate data for node Y=f(Xinv)+f(Xninv)+error
# input: sample size, variable x, error term (T/F), mean & sd of x, function types to use
# output: matrix containing values for (x,y)
nodedatagen = function(n = NULL, X_linear = matrix(numeric(0),0,0), X_nonlinear = matrix(numeric(0),0,0),
                       nodeN.linear = NULL, nodeN.nonlinear = NULL, se = 1, ninv.method = "gp",
                       inv.method = "random", data_fx="inv", error.distr = 'normal', edge_fx_list = NULL){

  X_linear = as.matrix(X_linear); X_nonlinear = as.matrix(X_nonlinear)
  nodeN.linear = ncol(X_linear)
  nodeN.nonlinear = ncol(X_nonlinear)
  # meant for linear edges
  if(nodeN.linear > 0) {
    Y_linear = sapply(1:nodeN.linear, function(i) pairdatagen(n, x = X_linear[, i], error = F, ninv.method = ninv.method,
                                                       inv.method = inv.method, se = se, data_fx=data_fx, error.distr = error.distr)$y)
  }else {
    Y_linear = matrix(0, n, 1)
  }
  # meant for non-linear edges
  if(length(edge_fx_list)==0){ # no pre-specified edge fx
    if(data_fx=="inv") edge_fx_list_nonlinear <- rep(inv.method, nodeN.nonlinear)
    if(data_fx=="ninv") edge_fx_list_nonlinear <- rep(ninv.method, nodeN.nonlinear)
    if(data_fx=="nonlinear") edge_fx_list_nonlinear <- rep("random", nodeN.nonlinear)
  }else{ # use pre-specified functions
    edge_fx_list_nonlinear <- edge_fx_list
  }

  if(nodeN.nonlinear > 0) {
    Y_nonlinear = sapply(1:nodeN.nonlinear, function(i) pairdatagen(n, x = X_nonlinear[, i], error = F, ninv.method = edge_fx_list_nonlinear[i],
                                                         inv.method = edge_fx_list_nonlinear[i], se = se, data_fx=data_fx, error.distr = error.distr)$y)
  }else{
    Y_nonlinear = matrix(0, n, 1)
  }

  se_val <- runif(1, min=se, max=se+0.25)

  if(is.null(error.distr)){
    e = rnorm(n, 0, se_val)
  }else if(error.distr == "gaussian"){
    e = rnorm(n, 0, se_val)
  }else if(error.distr == "laplace"){
    e = VGAM::rlaplace(n=n, location = 0, scale = se_val)
  }else if(error.distr == "gumbel"){
    e = VGAM::rgumbel(n=n, location = 0, scale = se_val)
  }else if(error.distr == "t"){
    e = extraDistr::rlst(n=n, df=5, mu=0, sigma = se_val)
  }

  y = rowSums(Y_linear) + rowSums(Y_nonlinear) + e

  dat = NULL
  dat$X_linear = X_linear
  dat$X_nonlinear = X_nonlinear
  dat$y = scale(y)
  return(dat)
}


#' Generates data according to DAG structure
#'
#' Given a DAG structure, the function generates according to the observed parent-child relations
#' Users may specify function types
#' If no DAG is provided, then function randomly generates a DAG
#'
#' @param dag_amat adjacency matrix of DAG
#' @param nodeN number of nodes in (random) DAG
#' @param dagSparsityProb DAG sparsity parameter
#' @param nonlinear.perct percentage of edges generated by invertible, non-linear functions
#' @param nonlinear.arcs list of arcs to be generated by non-linear functions
#' @param data_fx options: linear, inv, ninv, nonlinear(inv and ninv)
#' @param ninv.method non-invertible function to use
#' @param inv.method invertible function to use
#' @param n number of samples
#' @param labels variable names
#' @param se standard deviation value
#' @param edge_fx_amat adjacency matrix containing pre-specified functions for each edge
#' @param error.distr statistical distribution for error term
#'
#' @return DAG object of class bnlearn
#' @examples
#' amat <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
#' DAGdatagen(n=1000, dag_amat=amat, nodeN=rnow(amat))
#' @export

DAGdatagen = function(n = NULL, dag_amat = NULL, nodeN = nrow(dag_amat), labels = colnames(dag_amat), dagSparsityProb = NULL,
                      data_fx = "inv", nonlinear.perct = 1, nonlinear.arcs = NULL, se = 1,
                      ninv.method = "gp", inv.method = "random", error.distr='gaussian', edge_fx_amat=NULL) {
  # assign variable names
  if(is.null(labels)) labels = paste("X", 1:nodeN, sep = "")

  # create bnlearn obj from adjmat
  if(is.null(dag_amat)){
    if(is.null(dagSparsityProb)) dagSparsityProb = 2/(nodeN-1)
    dag.bn = bnlearn::random.graph(as.character(1:nodeN), prob = dagSparsityProb)
    bnlearn::nodes(dag.bn) = labels
  }else{
    colnames(dag_amat) <- labels
    dag.bn = bnlearn::empty.graph(colnames(dag_amat))
    bnlearn::amat(dag.bn) = dag_amat
    nodeN = length(dag.bn$nodes)
    labels = names(dag.bn$nodes)
  }


  arcs.ind = dag.bn$arcs
  n.arcs = nrow(arcs.ind)

  # non-linear edge percentage, default 100%
  n.arcs.nonlinear = round(n.arcs *  nonlinear.perct * (data_fx!="linear"))
  # declaring which are non-linear edges
  if(is.null(nonlinear.arcs)){
    nonlinear.arcs = matrix(arcs.ind[sample(1:n.arcs, n.arcs.nonlinear),], ncol= 2, byrow = F)
  }else{
    nonlinear.arcs = matrix(nonlinear.arcs, ncol = 2, byrow = F)
  }

  data = scale(sapply(1:nodeN, function(i)  rnorm(n, 0, 1)))
  colnames(data) = labels

  nodes.order = names(igraph::topo_sort(bnlearn::as.igraph(dag.bn)))

  for(i in nodes.order) {
    if(i %in% arcs.ind[, 2]) {
      if(i %in% nonlinear.arcs[, 2]) {
        pa.nonlinear= labels %in% nonlinear.arcs[nonlinear.arcs[, 2] %in% i, 1]
        pa.linear = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1][!arcs.ind[arcs.ind[, 2] %in% i, 1]%in%labels[pa.nonlinear]]
        if(!is.null(edge_fx_amat)){
          edge_fx_list <- edge_fx_amat[pa.nonlinear,i]
        }else{
          edge_fx_list <- NULL
        }
        data[, labels %in% i] = nodedatagen(n = n, X_linear = data[, pa.linear], X_nonlinear = data[, pa.nonlinear], se = se,
                                            ninv.method = ninv.method, inv.method = inv.method, data_fx=data_fx, error.distr=error.distr, edge_fx_list = edge_fx_list)$y
      }else{
        pa.linear = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1]
        data[, labels %in% i] = nodedatagen(n = n, X_linear = data[, pa.linear], se = se,
                                            ninv.method = ninv.method, inv.method = inv.method, data_fx="linear", error.distr=error.distr)$y
      }
    }
  }

  return(list(DAGdata = scale(data), nonlinear.arcs= nonlinear.arcs, dag.amat = bnlearn::amat(dag.bn)))
}


