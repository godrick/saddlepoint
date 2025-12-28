
test_that("linearlyMapped multinomial -> binomial: recovers exact MLE; using sparse A", {

  set.seed(2)
  B <- 12
  n_trials <- 8
  theta_true <- 0.35

  # Observed data: Y_b ~ Binomial(n_trials, theta_true)
  y <- rbinom(B, size = n_trials, prob = theta_true)
  y <- y[y!=0]
  B <- length(y)
  theta_hat_exact <- mean(y) / n_trials

  # a 3-category multinomial with probabilities:
  #   (theta/2, theta/2, 1-theta)
  prob_fun <- function(theta) c(theta[1] / 2, theta[1] / 2, 1 - theta[1])


  cgf_X <- MultinomialModelCGF(
    n = adaptor(fixed_param = n_trials),
    prob_vec = prob_fun,
    iidReps = 1
  )

  # Map X -> Y = X1 + X2 using a sparse matrix A
  A <- Matrix::Matrix(c(1, 1, 0), nrow = 1, sparse = TRUE)
  cgf_Y_block <- linearlyMappedCGF(cgf = cgf_X, matrix_A = A)

  # Replicate across B observations
  cgf_Y <- iidReplicatesCGF(cgf = cgf_Y_block, iidReps = B, block_size = 1)

  start_theta <- c(0.5)
  lb_theta <- c(1e-6)
  ub_theta <- c(1 - 1e-6)


  res_two <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = cgf_Y,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    method = "two_step"
  )

  res_con <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = cgf_Y,
    starting.theta = start_theta,
    lb.theta = lb_theta,
    ub.theta = ub_theta,
    std.error = TRUE,
    discrepancy = TRUE,
    method = "constrained"
  )

  expect_equal(as.numeric(res_two$MLEs.theta), theta_hat_exact, tolerance = 1e-3)
  expect_equal(as.numeric(res_con$MLEs.theta), theta_hat_exact, tolerance = 1e-3)

  # Correction term should not affect the 1-parameter binomial MLE; discrepancy should be ~ 0.
  expect_true(abs(as.numeric(res_two$discrepancy)) < 1e-6)
  expect_true(abs(as.numeric(res_con$discrepancy)) < 1e-6)

  # Both methods should agree
  expect_equal(res_con$MLEs.theta, res_two$MLEs.theta, tolerance = 1e-3)
})
