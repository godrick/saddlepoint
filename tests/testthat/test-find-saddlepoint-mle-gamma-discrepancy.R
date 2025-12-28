

.gamma_mle_exact <- function(x,
                             start = c(1, 1),
                             lower = c(1e-6, 1e-6),
                             upper = c(Inf, Inf)) {
  nll <- function(theta) -sum(dgamma(x, shape = theta[1], rate = theta[2], log = TRUE))
  nlminb(start = start, objective = nll, lower = lower, upper = upper)$par
}


test_that("Gamma: discrepancy is finite and (usually) moves SPA MLE toward exact MLE", {

  set.seed(303)
  n <- 10
  alpha_true <- 3.2
  beta_true <- 0.9
  x <- rgamma(n, shape = alpha_true, rate = beta_true)

  theta_exact <- .gamma_mle_exact(x)

  start_theta <- c(2, 1)
  lb_theta <- c(1e-6, 1e-6)
  ub_theta <- c(Inf, Inf)


  # Two-step fit
  res_two <- find.saddlepoint.MLE(
    observed.data = x,
    cgf = GammaCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    method = "two_step"
  )

  theta_spa <- res_two$MLEs.theta
  discr <- res_two$discrepancy

  expect_true(all(is.finite(theta_spa)))
  expect_true(all(is.finite(discr)))
  expect_true(max(abs(discr)) > 1e-10)  # typically non-zero for 2-parameter gamma

  # Discrepancy is intended to approximate (theta_exact - theta_spa).
  # so we check:
  #   - it points in the right direction, and
  #   - applying it as theta_spa + discrepancy reduces the error.
  delta_exact <- theta_exact - theta_spa
  err_before <- sqrt(sum(delta_exact^2))
  err_after  <- sqrt(sum((as.numeric(theta_exact) - (theta_spa + discr))^2))

  expect_true(sum(delta_exact * discr) > 0)
  expect_true(err_after <= err_before + 1e-8)

  # And it should be reasonably close relative to the true (exact-spa) shift.
  if (err_before > 1e-6) {
    rel_err <- sqrt(sum((delta_exact - discr)^2)) / err_before
    expect_true(rel_err < 0.5)
  } else {
    expect_true(max(abs(discr)) < 1e-6)
  }

  # Constrained fit should broadly agree with two-step for this small problem.
  res_con <- find.saddlepoint.MLE(
    observed.data = x,
    cgf = GammaCGF,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = FALSE,
    discrepancy = FALSE,
    method = "constrained"
  )

  # The two optimisation strategies should agree closely.
  expect_equal(res_con$MLEs.theta, res_two$MLEs.theta, tolerance = 1e-2)
})
