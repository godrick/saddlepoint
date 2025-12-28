
test_that("Poisson: SPA MLE equals exact MLE (mean); discrepancy ~ 0", {

  set.seed(202)
  n <- 12
  lambda_true <- 2.4
  y <- rpois(n, lambda = lambda_true)
  y <- y[y!=0]
  n <- length(y)

  lambda_hat_exact <- mean(y)

  start_theta <- c(1)
  lb_theta <- c(1e-8)
  ub_theta <- c(Inf)


  res_con <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = PoissonCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    zeroth.order = FALSE
  )

  expect_equal(res_con$MLEs.theta[1], lambda_hat_exact, tolerance = 1e-3)
  expect_true(abs(res_con$discrepancy[1]) < 1e-8)

  # Two-step should agree and uses nlminb
  res_two <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = PoissonCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    zeroth.order = FALSE,
    method = "two_step"
  )

  expect_equal(res_two$MLEs.theta[1], lambda_hat_exact, tolerance = 1e-3)
  expect_true(abs(res_two$discrepancy[1]) < 1e-8)
  expect_equal(res_two$MLEs.theta, res_con$MLEs.theta, tolerance = 1e-3)

  # Zeroth-order should give the same lambda MLE in Poisson (correction term does not depend on lambda).
  res_zero <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = PoissonCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = FALSE,
    discrepancy = FALSE,
    opts.user = opts,
    zeroth.order = TRUE,
    method = "two_step"
  )
  expect_equal(res_zero$MLEs.theta[1], lambda_hat_exact, tolerance = 1e-7)
})
