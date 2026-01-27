

test_that("NormalCGF: basic univariate identities", {
  mu <- 2.0
  sigma <- 1.5
  theta <- c(mu, sigma)

  tvec <- c(0.10, -0.20, 0.00, 0.30)

  # K(t) = sum_i mu t_i + 0.5 sigma^2 t_i^2
  K_ref <- sum(mu * tvec + 0.5 * sigma^2 * tvec^2)
  expect_equal(NormalCGF$K(tvec, theta), K_ref, tolerance = 1e-12)

  # K1(t) = mu + sigma^2 t (elementwise)
  K1_ref <- mu + sigma^2 * tvec
  expect_equal(NormalCGF$K1(tvec, theta), K1_ref, tolerance = 1e-12)

  # K2(t) is diagonal with sigma^2
  K2 <- NormalCGF$K2(tvec, theta)
  K2m <- as.matrix(K2)
  expect_equal(diag(K2m), rep(sigma^2, length(tvec)), tolerance = 1e-12)
  expect_equal(sum(abs(K2m - diag(diag(K2m)))), 0, tolerance = 1e-12)

  # tilting_exponent(t) = K(t) - t^T K1(t) = -0.5 sigma^2 sum(t^2)
  te_ref <- -0.5 * sigma^2 * sum(tvec^2)
  te_val <- NormalCGF$.private_api$tilting_exponent(tvec, theta)
  expect_equal(as.numeric(te_val), te_ref, tolerance = 1e-12)
})


test_that("NormalCGF: analytic saddlepoint t-hat matches closed form", {
  set.seed(1)
  y <- rnorm(7, mean = 1.3, sd = 0.9)
  mu <- 1.3
  sigma <- 0.9
  theta <- c(mu, sigma)

  t_hat <- NormalCGF$analytic_tvec_hat(y, theta)
  expect_equal(t_hat, (y - mu) / sigma^2, tolerance = 1e-12)
})


# -------------------------------------------------------------------------
# NormalModelCGF
# -------------------------------------------------------------------------

test_that("NormalModelCGF: adaptor composition agrees with NormalCGF", {
  # Parameterise mu and sigma as functions of theta
  mu_fn <- function(theta) theta[1]
  sigma_fn <- function(theta) exp(theta[2])

  cgf <- NormalModelCGF(mu = mu_fn, sigma = sigma_fn, iidReps = "any")

  theta <- c(0.4, log(1.7))
  mu <- mu_fn(theta)
  sigma <- sigma_fn(theta)

  tvec <- c(-0.3, 0.0, 0.25)

  expect_equal(
    cgf$K(tvec, theta),
    NormalCGF$K(tvec, c(mu, sigma)),
    tolerance = 1e-12
  )
  expect_equal(
    cgf$K1(tvec, theta),
    NormalCGF$K1(tvec, c(mu, sigma)),
    tolerance = 1e-12
  )
})


# -------------------------------------------------------------------------
# MultivariateNormalModelCGF
# -------------------------------------------------------------------------

test_that("MultivariateNormalModelCGF: replication semantics and derivatives", {
  d <- 2L
  B <- 3L

  # mu(theta) estimates a 2-vector; Sigma is fixed .
  mu_fn <- function(theta) theta[1:d]
  Sigma_fixed <- diag(c(1.0, 2.0))
  sigma_fn <- function(theta) Sigma_fixed

  mvn <- MultivariateNormalModelCGF(mu = mu_fn, sigma = sigma_fn, iidReps = B)

  theta <- c(0.2, -0.4)
  mu <- mu_fn(theta)
  Sigma <- Sigma_fixed

  tvec <- c(0.1, -0.2,
            0.0,  0.3,
           -0.1,  0.2)  # length = d*B

  tmat <- matrix(tvec, nrow = d)

  # K should be sum over blocks of mu^T t + 0.5 t^T Sigma t
  block_vals <- as.vector(crossprod(mu, tmat)) + 0.5 * colSums(tmat * (Sigma %*% tmat))
  K_ref <- sum(block_vals)
  expect_equal(as.numeric(mvn$K(tvec, theta)), K_ref, tolerance = 1e-12)

  # K1 is concatenation of mu + Sigma t for each block
  K1_ref <- as.vector(mu + Sigma %*% tmat)
  expect_equal(as.numeric(mvn$K1(tvec, theta)), K1_ref, tolerance = 1e-12)

  # K2 is block diagonal with Sigma on each diagonal block
  K2 <- as.matrix(mvn$K2(tvec, theta))
  expect_equal(K2[1:d, 1:d], Sigma, tolerance = 1e-12)
  expect_equal(K2[(d+1):(2*d), (d+1):(2*d)], Sigma, tolerance = 1e-12)
  expect_equal(K2[(2*d+1):(3*d), (2*d+1):(3*d)], Sigma, tolerance = 1e-12)
  expect_equal(max(abs(K2[1:d, (d+1):(3*d)])), 0, tolerance = 1e-12)

  # Analytic t-hat acts blockwise: Sigma^{-1} (y - mu)
  y <- rnorm(d * B)
  ymat <- matrix(y, nrow = d)
  t_hat_ref <- as.vector(solve(Sigma, ymat - mu))
  t_hat_val <- mvn$analytic_tvec_hat(y, theta)
  expect_equal(as.numeric(t_hat_val), t_hat_ref, tolerance = 1e-10)

  # logdetK2 should sum over blocks
  ld_ref <- B * as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  ld_val <- mvn$logdetK2(tvec, theta)
  expect_equal(as.numeric(ld_val), ld_ref, tolerance = 1e-10)

  # K2_solve should solve blockwise
  rhs <- rnorm(d * B)
  rhs_mat <- matrix(rhs, nrow = d)
  sol_ref <- as.vector(solve(Sigma, rhs_mat))
  sol_val <- mvn$K2_solve(tvec, theta, rhs)
  expect_equal(as.numeric(sol_val), sol_ref, tolerance = 1e-10)
})


test_that("MultivariateNormalModelCGF: errors for bad lengths and non-symmetric Sigma", {
  d <- 2L

  # Non-symmetric covariance (user error)
  mu_fn <- function(theta) theta[1:d]
  sigma_bad <- function(theta) matrix(c(1, 2,
                                        0, 1), nrow = 2, byrow = TRUE)

  mvn_bad <- MultivariateNormalModelCGF(mu = mu_fn, sigma = sigma_bad, iidReps = 1)

  theta <- c(0, 0)
  expect_error(
    mvn_bad$K(c(0.1, -0.2), theta),
    "Sigma must be symmetric"
  )

  # Bad tvec length for iidReps=1 (needs length == d)
  mvn_ok <- MultivariateNormalModelCGF(mu = mu_fn, sigma = function(th) diag(d), iidReps = 1)
  ## expect_error
  # mvn_ok$K(c(0.1, -0.2, 0.3), theta)

  # Bad tvec length for iidReps="any" (must be multiple of d)
  mvn_any <- MultivariateNormalModelCGF(mu = mu_fn, sigma = function(th) diag(d), iidReps = "any")
  expect_error(mvn_any$K(c(0.1, -0.2, 0.3), theta), "multiple")
})
















test_that("find.saddlepoint.MLE: NormalCGF (univariate) matches MVN(d=1) path", {

  set.seed(123)
  n <- 40
  mu_true <- 1.2
  sigma_true <- 0.8
  y <- rnorm(n, mean = mu_true, sd = sigma_true)

  mu_hat_exact <- mean(y)
  sigma_hat_exact <- sqrt(mean((y - mu_hat_exact)^2))

  start_theta <- c(mu_hat_exact + 0.25, sigma_hat_exact + 0.15)


  #  Univariate path
  res_uni <- find.saddlepoint.MLE(
    observed.data  = y,
    cgf            = NormalCGF,
    starting.theta = start_theta,
    lb.theta       = c(-Inf, 1e-3),
    std.error      = TRUE,
    discrepancy    = TRUE
  )

  #  Multivariate d=1 replicate path
  cgf_mv1 <- MultivariateNormalModelCGF(
    mu    = function(th) th[1],
    sigma = function(th) matrix(th[2]^2, nrow = 1),
    iidReps = n
  )

  res_mv1 <- find.saddlepoint.MLE(
    observed.data  = y,
    cgf            = cgf_mv1,
    starting.theta = start_theta,
    lb.theta       = c(-Inf, 1e-6),
    std.error      = TRUE,
    discrepancy    = TRUE,
    zeroth.order   = FALSE
  )

  # Both should land at the same MLEs (up to tolerance)
  expect_equal(res_uni$MLEs.theta, res_mv1$MLEs.theta, tolerance = 1e-6)

  # Both should match the true MLEs
  expect_equal(res_uni$MLEs.theta[1], mu_hat_exact, tolerance = 1e-4)
  expect_equal(res_uni$MLEs.theta[2], sigma_hat_exact, tolerance = 1e-4)

  expect_equal(res_mv1$MLEs.theta[1], mu_hat_exact, tolerance = 1e-4)
  expect_equal(res_mv1$MLEs.theta[2], sigma_hat_exact, tolerance = 1e-4)

  expect_true(all(is.finite(res_mv1$inverse.hessian)))
  expect_true(all(is.finite(res_uni$inverse.hessian)))


  if (all(is.finite(res_mv1$inverse.hessian)) && all(is.finite(res_uni$inverse.hessian))) {
    expect_equal(
      as.numeric(res_uni$inverse.hessian),
      as.numeric(res_mv1$inverse.hessian),
      tolerance = 1e-5
    )
  }
})


test_that("find.saddlepoint.MLE: non-symmetric Sigma triggers an error", {

  d <- 2L

  mu_fn <- function(theta) theta[1:d]
  sigma_bad <- function(theta) matrix(c(1, 2,
                                        0, 1), nrow = 2, byrow = TRUE) # NOT symmetric

  mvn_bad <- MultivariateNormalModelCGF(mu = mu_fn, sigma = sigma_bad, iidReps = 1)

  set.seed(1)
  y <- rnorm(d)
  start_theta <- c(0.7, 0.5)

  # We only need it to *error*, not to optimize
  expect_error(
    find.saddlepoint.MLE(
      observed.data  = y,
      cgf            = mvn_bad,
      starting.theta = start_theta,
      std.error      = TRUE,
      discrepancy    = TRUE,
      opts.user      = list(maxeval = 5, print_level = 0)
    )
  )
})
