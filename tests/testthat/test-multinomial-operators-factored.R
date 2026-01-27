
make_Q_factored <- function(d, seed = 1) {
  set.seed(seed)
  M <- matrix(rnorm(d * d), d, d)
  Q <- crossprod(M) + diag(0.1, d)
  cholQ <- chol(Q)
  dd <- diag(cholQ)^2
  A <- t(cholQ) %*% diag(1 / diag(cholQ))
  list(Q = Q, A = A, d = dd)
}


test_that("K4operatorAABB_factored matches explicit K4operatorAABB for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)
  tvec <- c(0.05, -0.08, 0.02)

  fd <- make_Q_factored(3, seed = 10)

  K4_fact <- cgf$.private_api$K4operatorAABB_factored

  got_fact <- K4_fact(tvec, param, fd$A, fd$d, fd$A, fd$d)
  got_Q    <- cgf$K4operatorAABB(tvec, param, fd$Q, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})


test_that("K3K3operatorAABBCC_factored matches explicit K3K3operatorAABBCC for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 8
  odds <- c(0.2, 0.3, 0.5)
  param <- c(N, odds)
  tvec <- c(0.02, -0.03, 0.01)

  fd <- make_Q_factored(3, seed = 20)

  K3K3_fact <- cgf$.private_api$K3K3operatorAABBCC_factored

  got_fact <- K3K3_fact(tvec, param, fd$A, fd$d, fd$A, fd$d, fd$A, fd$d)
  got_Q    <- cgf$K3K3operatorAABBCC(tvec, param, fd$Q, fd$Q, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})


test_that("K3K3operatorABCABC_factored matches explicit K3K3operatorABCABC for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 9
  odds <- c(3, 1, 6)
  param <- c(N, odds)
  tvec <- c(0.04, -0.01, -0.02)

  fd <- make_Q_factored(3, seed = 30)

  K3K3_fact <- cgf$.private_api$K3K3operatorABCABC_factored

  got_fact <- K3K3_fact(tvec, param, fd$A, fd$d, fd$A, fd$d, fd$A, fd$d)
  got_Q    <- cgf$K3K3operatorABCABC(tvec, param, fd$Q, fd$Q, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})
