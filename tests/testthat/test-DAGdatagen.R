

test_that("data dimensions are correct", {
  data <- DAGdatagen(n = 200, nodeN = 5, data_fx='nonlinear', ninv.method="sigmoid", se=1)$DAGdata
  expect_equal(dim(data), c(200, 5))
})


test_that("data dimensions are correct, given adjmat", {
  adjmat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0),
                  c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
  data <- DAGdatagen(n = 100, dag_amat = adjmat, data_fx='nonlinear', ninv.method="sigmoid", se=1)$DAGdata
  expect_equal(dim(data), c(100, 6))
})
