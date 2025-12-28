
# -------------------------------------------------------------------------
# Univariate Normal: saddlepoint MLE should match exact MLE
# -------------------------------------------------------------------------

test_that("find.saddlepoint.MLE: univariate Normal matches exact MLE (standard)", {
  set.seed(123)

  n <- 40
  mu_true <- 1.2
  sigma_true <- 0.8
  y <- rnorm(n, mean = mu_true, sd = sigma_true)

  mu_hat_exact <- mean(y)
  sigma_hat_exact <- sqrt(mean((y - mu_hat_exact)^2))

  # Saddlepoint MLE
  start_theta <- c(mu_hat_exact + 0.25, sigma_hat_exact + 0.15)

  res <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = NormalCGF,
    starting.theta = start_theta,
    lb.theta = c(-Inf, 0.001),
    lb.tvec = rep(-Inf, n),
    ub.tvec = rep( Inf, n),
    std.error = TRUE,
    discrepancy = TRUE
  )

  expect_true(all(is.finite(res$MLEs.theta)))
  expect_equal(res$MLEs.theta[1], mu_hat_exact, tolerance = 1e-3)
  expect_equal(res$MLEs.theta[2], sigma_hat_exact, tolerance = 1e-3)

  expect_true(all(is.finite(res$std.error)))

  expect_true(all(is.finite(res$discrepancy)))
})


# test_that("find.saddlepoint.MLE: univariate Normal runs under zeroth-order objective", {
#   set.seed(123)
#
#   n <- 30
#   y <- rnorm(n, mean = 0.3, sd = 1.1)
#
#   res0 <- find.saddlepoint.MLE(
#     observed.data = y,
#     cgf = NormalCGF,
#     starting.theta = c(0, 1),
#     lb.theta = c(-Inf, 0.0001),
#     ub.theta = c( Inf, 1e6),
#     std.error = FALSE,
#     discrepancy = TRUE,
#     zeroth.order = TRUE
#   )
#
#   expect_true(all(is.finite(res0$MLEs.theta)))
#   expect_true(all(is.finite(res0$discrepancy)))
# })


# -------------------------------------------------------------------------
# Multivariate Normal: estimate mu with fixed Sigma
# -------------------------------------------------------------------------

test_that("find.saddlepoint.MLE: multivariate Normal (mu only, Sigma fixed)", {
  set.seed(42)

  d <- 2
  n <- 25

  mu_true <- c(-0.5, 0.7)
  Sigma <- matrix(c(1.0, 0.3,
                    0.3, 2.0), nrow = d, byrow = TRUE)

  # Simulate n i.i.d. draws
  L <- chol(Sigma)
  Y <- matrix(rnorm(d * n), nrow = d)
  Y <- mu_true + L %*% Y  # d x n

  # Observed data accepted as matrix: columns = replicate blocks
  y_obs <- Y

  # Exact MLE for mu with known Sigma is the sample mean
  mu_hat_exact <- rowMeans(Y)

  # Model: theta is mu (length d), Sigma fixed
  mvn_cgf <- MultivariateNormalModelCGF(
    mu = adaptor(indices = 1:d),
    sigma = function(theta) Sigma,
    iidReps = "any"
  )

  res <- find.saddlepoint.MLE(
    observed.data = y_obs,
    cgf = mvn_cgf,
    starting.theta = c(0, 0.0005),
    lb.theta = c(-Inf, 0.0001),
    std.error = TRUE,
    discrepancy = TRUE
  )

  expect_equal(res$MLEs.theta, mu_hat_exact, tolerance = 1e-3)
  expect_true(all(is.finite(res$std.error)))
  expect_true(all(is.finite(res$discrepancy)))
})
