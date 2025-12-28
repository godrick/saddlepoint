
test_that("RSS Poisson+Bernoulli: exact likelihood implies p_hat = mean(y)/lambda; discrepancy ~ 0", {

  set.seed(303)
  B <- 10
  lambda_N <- 4
  p_true <- 0.30

  # Under Poisson thinning, Y ~ Poisson(lambda_N * p)
  y <- rpois(B, lambda = lambda_N * p_true)
  y <- y[y!=0]
  B <- length(y)
  p_hat_exact <- mean(y) / lambda_N
  p_hat_exact <- max(min(p_hat_exact, 1), 0)

  # Build RSS CGF: N ~ Poisson(lambda_N) and X ~ Bernoulli(p)
  count_cgf <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda_N), iidReps = 1)
  summand_cgf <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(indices = 1),
    iidReps = 1
  )
  rss_cgf <- randomlyStoppedSumCGF(count_cgf = count_cgf, summand_cgf = summand_cgf, iidReps = B)

  start_theta <- c(0.5)
  lb_theta <- c(1e-6)
  ub_theta <- c(1 - 1e-6)
  opts <- list(maxeval = 300, xtol_rel = 1e-10, print_level = 0)

  # two_step
  res_two <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = rss_cgf,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts,
    method = "two_step"
  )

  expect_equal(res_two$MLEs.theta, p_hat_exact, tolerance = 2e-3)
  expect_true(is.finite(res_two$discrepancy))
  expect_lt(abs(res_two$discrepancy), 5e-4)

  # constrained
  res_con <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = rss_cgf,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts,
    method = "constrained"
  )

  expect_equal(as.numeric(res_con$MLEs.theta), p_hat_exact, tolerance = 2e-3)
  expect_lt(abs(as.numeric(res_con$discrepancy)), 5e-4)

  # Matrix observed.data should be flattened by columns without changing the result.
  y_mat <- matrix(y, nrow = 1)
  expect_message(
    res_mat <- find.saddlepoint.MLE(
      observed.data = y_mat,
      cgf = rss_cgf,
      starting.theta = start_theta,
      lb.theta = lb_theta,
      ub.theta = ub_theta,
      std.error = FALSE,
      discrepancy = FALSE,
      opts.user = opts,
      method = "two_step"
    ),
    "Treating columns"
  )
  expect_equal(as.numeric(res_mat$MLEs.theta), as.numeric(res_two$MLEs.theta), tolerance = 1e-10)
})
