
test_that("randomlyStoppedSumCGF: ineq_constraint composes count-domain and summand-domain constraints", {
  # Count: Geometric(prob) has a scalar domain constraint in its CGF argument.
  count_cgf <- GeometricModelCGF(prob = adaptor(fixed_param = 0.4), iidReps = 1)

  # Summand: Gamma(shape, rate) has constraint t < rate.
  summand_cgf <- GammaModelCGF(shape = adaptor(fixed_param = 2.0),
                               rate  = adaptor(fixed_param = 3.0),
                               iidReps = 1)

  cgf <- randomlyStoppedSumCGF(count_cgf = count_cgf,
                               summand_cgf = summand_cgf,
                               iidReps = "any",
                               block_size = 1)

  param <- c(1.0)  # dummy param vector

  t <- 0.1
  Kx <- summand_cgf$K(t, param)

  ref <- c(count_cgf$ineq_constraint(Kx, param),
           summand_cgf$ineq_constraint(t, param))

  got <- cgf$ineq_constraint(t, param)

  expect_equal(as.numeric(got), as.numeric(ref), tolerance = 1e-12)

  # At a clearly feasible t, all constraints should be strictly negative
  expect_true(all(got < 0))
})
