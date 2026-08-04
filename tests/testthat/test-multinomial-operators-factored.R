make_Q_factored <- function(d, rank = max(1L, d - 1L), seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(d * rank), nrow = d, ncol = rank)
  dd <- runif(rank, min = 0.4, max = 1.2)
  Q <- A %*% diag(dd, nrow = rank) %*% t(A)
  list(Q = Q, A = A, d = dd)
}

expect_multinomial_factored_matches_dense <- function(cgf, param, tvec, seed) {
  fd <- make_Q_factored(length(tvec), rank = 2, seed = seed)

  expect_equal(
    as.numeric(cgf$.private_api$K4operatorAABB_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K4operatorAABB(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorAABBCC_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K3K3operatorAABBCC(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorABCABC_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K3K3operatorABCABC(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
}

test_that("multinomial factored contractions match dense contractions", {
  cgf <- MultinomialFamilyCGF$new()
  expect_true(saddlepoint:::.factored_delegate_is_safe(
    cgf$.private_api$K3K3operatorAABBCC_factored
  ))
  expect_true(saddlepoint:::.factored_delegate_is_safe(
    cgf$.private_api$K3K3operatorABCABC_factored
  ))

  expect_multinomial_factored_matches_dense(
    cgf,
    param = c(10, 2, 3, 5),
    tvec = c(0.05, -0.08, 0.02),
    seed = 10
  )
})

test_that("thin multinomial ABCABC preserves values and two AD orders", {
  set.seed(135)
  n <- 20L
  rank <- 2L
  tvec <- seq(-0.08, 0.09, length.out = n)
  A0 <- matrix(rnorm(n * rank), n, rank) +
    outer(rep(1, n), c(1e4, -2e4))
  A1 <- matrix(rnorm(n * rank), n, rank)
  odds0 <- seq(0.6, 1.4, length.out = n)
  theta <- 0.15
  cgf <- MultinomialFamilyCGF$new()

  parameter <- function(x) {
    c(10 + exp(x[1]), odds0 * exp(x[1] * seq_len(n) / (20 * n)))
  }
  factor_A <- function(x) A0 + 0.02 * x[1] * A1
  factor_d <- function(x) c(exp(0.1 * x[1]), exp(-0.05 * x[1]))

  candidate <- RTMB::MakeTape(function(x) {
    cgf$.private_api$K3K3operatorABCABC_factored(
      tvec, parameter(x), factor_A(x), factor_d(x)
    )
  }, theta)
  reference <- RTMB::MakeTape(function(x) {
    A <- factor_A(x)
    d <- factor_d(x)
    param <- parameter(x)
    v <- param[-1] * exp(tvec)
    v <- v / sum(v)
    A <- A - outer(rep(1, n), as.vector(crossprod(v, A)))
    Q <- A %*% (d * t(A))
    cgf$K3K3operatorABCABC(tvec, param, Q)
  }, theta)

  expect_equal(candidate(theta), reference(theta), tolerance = 1e-8)
  expect_equal(
    candidate$jacobian(theta),
    reference$jacobian(theta),
    tolerance = 1e-7
  )
  expect_equal(
    candidate$jacfun()$jacobian(theta),
    reference$jacfun()$jacobian(theta),
    tolerance = 1e-6
  )
})

test_that("multinomial contractions remove null offsets and balance factors", {
  cgf <- MultinomialFamilyCGF$new()
  param <- c(7, 0.15, 0.25, 0.60)
  tvec <- c(0.03, -0.02, 0.01)
  u1 <- c(0.2, -0.4, 0.7)
  u2 <- c(-0.1, 0.5, 0.3)
  u3 <- c(0.6, 0.2, -0.2)
  u4 <- c(0.4, -0.3, 0.1)

  expect_equal(
    cgf$K3operator(tvec, param, u1 + 1e4, u2 - 2e4, u3 + 3e4),
    cgf$K3operator(tvec, param, u1, u2, u3),
    tolerance = 1e-9
  )
  expect_equal(
    cgf$K4operator(tvec, param, u1 + 1e4, u2 - 2e4, u3 + 3e4, u4 - 4e4),
    cgf$K4operator(tvec, param, u1, u2, u3, u4),
    tolerance = 1e-8
  )

  tvec2 <- c(0, 0)
  A <- matrix(c(1e250, 0), 2L, 1L)
  dvec <- 1e-320
  Q <- A %*% (dvec * t(A))
  parameter <- function(p) c(exp(p[1]), 1e-119, 1)
  candidate <- RTMB::MakeTape(function(p) {
    cgf$.private_api$K3K3operatorAABBCC_factored(
      tvec2, parameter(p), A, dvec
    )
  }, 0)
  reference <- RTMB::MakeTape(function(p) {
    cgf$K3K3operatorAABBCC(tvec2, parameter(p), Q)
  }, 0)
  actual <- c(
    candidate(0),
    candidate$jacobian(0),
    candidate$jacfun()$jacobian(0)
  )
  expected <- c(
    reference(0),
    reference$jacobian(0),
    reference$jacfun()$jacobian(0)
  )
  expect_true(all(is.finite(actual)))
  expect_equal(actual, expected, tolerance = 1e-12)
})
