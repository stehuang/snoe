
test_that("algorithm returns adjacency matrix", {
  data <- DAGdatagen(n = 200, nodeN = 5, data_fx='inv', ninv.method="random", se=1)$DAGdata
  adjmat <- nleo(data=data)
  expect_equal(dim(adjmat), c(5,5))
})


test_that("algorithm returns DAG", {
  data <- DAGdatagen(n = 200, nodeN = 10, data_fx='inv', ninv.method="random", se=1)$DAGdata
  adjmat <- nleo(data=data)
  dag_obj <- bnlearn::empty.graph(colnames(adjmat))
  amat(dag_obj) <- adjmat
  expect_equal(acyclic(dag_obj), T)
})
