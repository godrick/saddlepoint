test_that("PoissonCGF basic functions match definitions", {
  tvec <- c(-0.2, 0, 0.3, 1.1)
  lambda <- 2
  expect_equal(PoissonCGF$K(tvec, lambda), lambda * (exp(tvec) - 1))
  expect_equal(PoissonCGF$K1(tvec, lambda), lambda * exp(tvec))
  expect_equal(PoissonCGF$K2(tvec, lambda), lambda * exp(tvec))
  expect_equal(PoissonCGF$K3operator(tvec, lambda, rep(1, length(tvec)), rep(1, length(tvec)), rep(1, length(tvec))),
               sum(lambda * exp(tvec)))
  expect_equal(PoissonCGF$K4operator(tvec, lambda, rep(1, length(tvec)), rep(1, length(tvec)), rep(1, length(tvec)), rep(1, length(tvec))),
               sum(lambda * exp(tvec)))
  # analytic t-hat
  x <- c(1.5, 2.0, 3.5, 4.0)
  expect_equal(PoissonCGF$analytic_tvec_hat(x, lambda), log(x / lambda))
})

test_that("PoissonCGF broadcasting with iidReps='any' works", {
  tvec <- c(-0.1, 0.2, 0.3, -0.4)
  lam <- c(2, 5)
  # Expect broadcasting to [2,5,2,5]
  expected_lam <- rep(lam, length.out = length(tvec))
  expect_equal(PoissonCGF$K1(tvec, lam), expected_lam * exp(tvec))
  expect_equal(PoissonCGF$K(tvec, lam), expected_lam * (exp(tvec) - 1))
})

test_that("PoissonModelCGF with lambda function works (iidReps='any')", {
  lambda_fn <- function(theta) c(theta[1], 3 * theta[1])
  K <- PoissonModelCGF(lambda = lambda_fn, iidReps = "any")
  theta0 <- 2
  tvec <- c(0.0, 0.1, -0.2, 0.3)
  lam <- lambda_fn(theta0)
  expected_lam <- rep(lam, length.out = length(tvec))
  expect_equal(K$K1(tvec, theta0), expected_lam * exp(tvec))
  # analytic t-hat
  x <- c(0.7, 2.5, 1.3, 3.0)
  expect_equal(K$analytic_tvec_hat(x, theta0), log(x / expected_lam))
})

test_that("PoissonModelCGF enforces iidReps when integer is provided", {
  K <- PoissonModelCGF(lambda = function(theta) c(3, 4), iidReps = 3)
  # valid: length(tvec) == 2 * 3
  t_ok <- rep(0, 6)
  expect_silent(K$K1(t_ok, 1))
  # invalid: wrong length
  t_bad <- rep(0, 5)
  expect_error(K$K1(t_bad, 1))
})

