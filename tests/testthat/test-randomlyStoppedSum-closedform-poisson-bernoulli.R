
test_that("randomlyStoppedSumCGF: Poisson count + Bernoulli(=Binomial(1,p)) summand reduces to Poisson(lambda*p)", {
  # Count: N ~ Poisson(lambda)
  count_cgf <- PoissonModelCGF(lambda = adaptor(indices = 2), iidReps = 1)

  # Summand: X ~ Bernoulli(p) as Binomial(n=1,p)
  summand_cgf <- BinomialModelCGF(n = adaptor(fixed_param = 1),
                                  p = adaptor(indices = 1),
                                  iidReps = 1)

  # Stopped sum: Y = sum_{i=1}^N X_i (scalar)
  cgf <- randomlyStoppedSumCGF(count_cgf = count_cgf,
                               summand_cgf = summand_cgf,
                               iidReps = "any",
                               block_size = 1)

  theta <- c(p = 0.3, lambda = 4.2)
  t     <- 0.13

  mu <- theta[["lambda"]] * theta[["p"]]  # mean for the implied Poisson

  K_ref  <- mu * (exp(t) - 1)
  K1_ref <- mu * exp(t)
  K2_ref <- mu * exp(t)
  K3_ref <- mu * exp(t)
  K4_ref <- mu * exp(t)

  expect_equal(as.numeric(cgf$K(t, theta)),  K_ref,  tolerance = 1e-12)
  expect_equal(as.numeric(cgf$K1(t, theta)), K1_ref, tolerance = 1e-10)
  expect_equal(as.numeric(cgf$K2(t, theta)), K2_ref, tolerance = 1e-10)

  set.seed(1)
  w <- rnorm(4)

  expect_equal(
    as.numeric(cgf$K2operator(t, theta, w[1], w[2])),
    K2_ref * w[1] * w[2],
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(cgf$K3operator(t, theta, w[1], w[2], w[3])),
    K3_ref * w[1] * w[2] * w[3],
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(cgf$K4operator(t, theta, w[1], w[2], w[3], w[4])),
    K4_ref * w[1] * w[2] * w[3] * w[4],
    tolerance = 1e-7
  )

  # K2operatorAK2AT sanity (A is 2x1)
  A <- matrix(rnorm(2), nrow = 2)
  ref_A <- A %*% matrix(K2_ref, 1, 1) %*% t(A)
  got_A <- cgf$K2operatorAK2AT(t, theta, A)
  expect_equal(as.matrix(got_A), as.matrix(ref_A), tolerance = 1e-8)

  # ineq_constraint should just be empty here (both Poisson & Binomial have full domain)
  expect_equal(length(cgf$ineq_constraint(t, theta)), 0L)
})
