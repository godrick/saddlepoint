
test_that("RSS: independent replicates vs shared-count multivariate model give different MLEs", {

  B <- 5
  lambda_N <- 4

  # A deliberately heterogeneous observation vector:
  # under an i.i.d. model this is fine; under a shared-count model it pushes p.
  y <- c(5, 4, 1, 5, 1)

  # Exact MLE under the i.i.d. model (Poisson thinning): p_hat = mean(y)/lambda_N
  p_hat_iid_exact <- mean(y) / lambda_N

  count_cgf <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda_N), iidReps = 1)
  summand_scalar <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(indices = 1),
    iidReps = 1
  )

  # Correct model: B independent replicates of scalar RSS
  cgf_iid <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_scalar,
    iidReps = B
  )

  # Wrong-but-plausible model: one B-dimensional RSS with a shared count N,
  # achieved by treating the *summand* as a B-vector (independent components) and
  # then summing over a common N.
  summand_vec <- iidReplicatesCGF(cgf = summand_scalar, iidReps = B, block_size = 1)
  cgf_sharedN <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_vec,
    iidReps = 1
  )

  opts <- list(maxeval = 300, xtol_rel = 1e-10, print_level = 0)

  fit_iid <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = cgf_iid,
    starting.theta = 0.2,
    lb.theta = 1e-6,
    ub.theta = 1 - 1e-6,
    method = "two_step",
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts
  )

  fit_shared <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = cgf_sharedN,
    starting.theta = 0.2,
    lb.theta = 1e-6,
    ub.theta = 1 - 1e-6,
    method = "two_step",
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts
  )

  # i.i.d. model: saddlepoint MLE should equal the exact MLE (Poisson thinning)
  expect_equal(fit_iid$MLEs.theta, p_hat_iid_exact, tolerance = 1e-4)

  # shared-count model should generally produce a *different* MLE
  expect_true(abs(as.numeric(fit_shared$MLEs.theta) - as.numeric(fit_iid$MLEs.theta)) > 0.02)
})
