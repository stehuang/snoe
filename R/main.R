

#' Nonlinear Causal Learning Algorithm by Likelihood Ratio Test
#'
#' Function for running the full algorithm
#'
#' @param data full dataset
#' @param skel_alpha sig. level for learning the initial skeleton
#' @param test_alpha sig. level for likelihood ratio test
#' @param covar_alpha sig. level for covariate testing
#' @param ll_test_approach variation of likelihood test to use; options are "2fold" and "5_5"
#' @param train_ratio percentage of data used for training dataset
#' @param standardized standardize data or not?
#'
#' @return adjacency matrix of learned DAG
#' @export
#'
#' @examples
#' nleo(data)
nleo <- function(data, skel_alpha=0.25, test_alpha=0.05, covar_alpha=0.00001,
                 ll_test_approach="5_5", train_ratio=0.5, standardized=T){
  if(standardized){
    data <- apply(data, 2, function(x) if(is.numeric(x)){scale(x, center=TRUE, scale=TRUE)} else x)
  }
  data <- as.data.frame(data)

  # 1. estimate skeleton
  init_amat <- estimate_cpdag(data = data, learning_alg = "PC", pc_params = list(alpha=test_alpha, alpha_dense=skel_alpha, test='cor'))
  print(init_amat)
  # 2. edge orientation on cpdag
  oriented_amat <- orient_edge(data = data, curr_amat=init_amat, alpha=test_alpha,
                                train_ratio = train_ratio, ll_test_approach = ll_test_approach)
  print(oriented_amat)
  # 3. edge deletion
  pruned_amat <- prune_dag(data = data, curr_amat=oriented_amat, alpha=covar_alpha)
  print(pruned_amat)
  # 4. extend PDAG to DAG
  final_amat <- finalize_dag(data = data, curr_amat=pruned_amat, alpha=test_alpha, train_ratio = train_ratio, ll_test_approach = ll_test_approach)
  print(final_amat)
  return(final_amat)
}
