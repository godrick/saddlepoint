

test_that("Normal: both methods recover the exact MLE; discrepancy ~ 0", {

  set.seed(101)
  n <- 20
  mu_true <- 5
  sigma_true <- 2
  y <- rnorm(n, mean = mu_true, sd = sigma_true)

  # Exact MLE for Normal(mu, sigma)
  mu_hat_exact <- mean(y)
  sigma_hat_exact <- sqrt(mean((y - mu_hat_exact)^2))

  start_theta <- c(0, 1)
  lb_theta <- c(-Inf, 1e-6)
  ub_theta <- c( Inf,  Inf)


  res_con <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = NormalCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    starting.tvec = rep(0, length(y)),
    lb.tvec = rep(-Inf, length(y)),
    ub.tvec = rep( Inf, length(y)),
    std.error = TRUE,
    discrepancy = TRUE,
    method = "constrained"
  )

  expect_equal(res_con$MLEs.theta[1], mu_hat_exact, tolerance = 1e-5)
  expect_equal(res_con$MLEs.theta[2], sigma_hat_exact, tolerance = 1e-5)
  expect_true(max(abs(res_con$discrepancy)) < 1e-6)

  # Analytic t-hat for Normal: t_i = (y_i - mu) / sigma^2
  t_hat_exact <- (y - res_con$MLEs.theta[1]) / (res_con$MLEs.theta[2]^2)
  expect_equal(res_con$MLEs.tvec, t_hat_exact, tolerance = 1e-6)

  res_two <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = NormalCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    starting.tvec = rep(0, length(y)),
    lb.tvec = rep(-Inf, length(y)),
    ub.tvec = rep( Inf, length(y)),
    std.error = TRUE,
    discrepancy = TRUE,
    method = "two_step",
    control.nlminb = list(eval.max = 200, iter.max = 150, rel.tol = 1e-10)
  )

  expect_equal(res_two$MLEs.theta, res_con$MLEs.theta, tolerance = 1e-5)
  expect_equal(res_two$MLEs.tvec,  res_con$MLEs.tvec,  tolerance = 1e-3)
  expect_true(isTRUE(res_two$optimizer == "nlminb"))
  expect_true(max(abs(res_two$discrepancy)) < 1e-6)
})
