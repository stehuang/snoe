% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_cpdag.R
\name{estimate_cpdag}
\alias{estimate_cpdag}
\title{STEP 1 OF ALGORITHM: Learning the initial CPDAG from data}
\usage{
estimate_cpdag(
  data,
  nodelabels = colnames(data),
  cluster = NULL,
  learning_alg = "PC",
  pc_params = list(alpha = 0.01, alpha_dense = 0.25, test = "cor"),
  ccdr_pc_alpha = 0.01,
  ccdr_params = list(pc_alpha = ccdr_pc_alpha, lambdas.length = 15, gamma = 2, max.iters
    = sparsebnUtils::default_max_iters(ncol(data)), ccdr_alpha =
    sparsebnUtils::default_alpha(), error.tol = 1e-04, verbose = FALSE),
  debug = FALSE
)
}
\arguments{
\item{data}{observed data}

\item{nodelabels}{vector of variable names}

\item{cluster}{use cluster for computing or not (bnlearn function argument)}

\item{learning_alg}{options are "PC" or "ccdr"}

\item{pc_params}{parameters to pass into PC algorithm; alpha is threshold for CI tests for v-structure detection, alpha_dense is for learning the skeleton}

\item{ccdr_pc_alpha}{estimates CPDAG from PC to determine size of DAG from ccdr}

\item{ccdr_params}{parameters to pass into ccdr algorithm}

\item{debug}{display debug process or not (bnlearn function argument)}
}
\value{
adjacency matrix of learned cpdag
}
\description{
Learns a dense CPDAG from data; skeleton is estimated under a relaxed threshold, and conditional indep relations are detected using a strigent threshold
}
\examples{
estimate_cpdag(data, standardized=T, learning_alg="PC", pc_params(alpha=0.05, alpha_dense=0.25, test='cor'))
}
