test_that("Matrix observed.data is flattened column-wise and treated as iid replicate blocks", {

  set.seed(505)
  d <- 2
  B <- 20

  mu_true <- c(1.0, -4.5)
  Sigma <- matrix(c(1.0, 0.3, 0.3, 2.0), nrow = 2, byrow = TRUE)
  U <- chol(Sigma)  # upper-triangular, t(U)%*%U = Sigma

  Z <- matrix(rnorm(d * B), nrow = d)
  Ymat <- matrix(mu_true, nrow = d, ncol = B) + t(U) %*% Z

  # Exact MLE of mu with known Sigma
  mu_hat_exact <- rowMeans(Ymat)

  cgf <- MultivariateNormalModelCGF(
    mu = adaptor(indices = 1:d),
    sigma = adaptor(fixed_param = Sigma),
    iidReps = B
  )

  start_theta <- c(0, 0)

  res_mat <- expect_message(
    find.saddlepoint.MLE(
      observed.data = Ymat,
      cgf = cgf,
      starting.theta = start_theta,
      std.error = TRUE,
      discrepancy = TRUE,
      method = "two_step"
    ),
    "Treating columns"
  )

  res_vec <- find.saddlepoint.MLE(
    observed.data = as.numeric(Ymat),
    cgf = cgf,
    starting.theta = start_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts,
    method = "two_step"
  )

  expect_equal(res_vec$MLEs.theta, mu_hat_exact, tolerance = 1e-6)

  # Normal is exact under SPA, so discrepancy should be ~0
  expect_true(max(abs(res_vec$discrepancy)) < 1e-17)
})
