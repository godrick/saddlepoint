# NOTE: MultinomialCGF is created via createMultinomialFamilyCGF(),
# which wraps a single-block MultinomialFamilyCGF in iidReplicatesCGF().
# These tests verify that that wrapper aggregates blocks correctly.



test_that("Multinomial iidReplicates wrapper aggregates blocks correctly", {
  # Single-block reference CGF (expects length(tvec) == d)
  base <- MultinomialFamilyCGF$new(iidReps = "any")

  # Aggregated CGF (wraps base in iidReplicatesCGF)
  agg  <- createMultinomialFamilyCGF(iidReps = "any", op_name = "MultinomialCGF_test")

  N <- 7
  odds <- c(2, 1, 4)
  param <- c(N, odds)

  d <- length(odds)
  B <- 4L

  t_block <- c(0.05, -0.03, 0.01)
  tvec <- rep(t_block, B)

  # ---- K: sum over blocks ----
  ref_K <- B * base$K(t_block, param)
  got_K <- agg$K(tvec, param)
  expect_equal(as.numeric(got_K), as.numeric(ref_K), tolerance = 1e-12)

  # ---- K1: concatenation ----
  ref_K1 <- rep(as.numeric(base$K1(t_block, param)), B)
  got_K1 <- as.numeric(agg$K1(tvec, param))
  expect_equal(got_K1, ref_K1, tolerance = 1e-12)

  # ---- K2: block diagonal ----
  K2_block <- as.matrix(base$K2(t_block, param))
  K2_ref <- matrix(0, nrow = d * B, ncol = d * B)
  for (b in seq_len(B)) {
    idx <- ((b - 1L) * d + 1L):(b * d)
    K2_ref[idx, idx] <- K2_block
  }
  got_K2 <- as.matrix(agg$K2(tvec, param))
  expect_equal(got_K2, K2_ref, tolerance = 1e-12)
})



test_that("Multinomial wrapper enforces fixed iidReps when requested", {
  base <- MultinomialFamilyCGF$new(iidReps = "any")

  N <- 7
  odds <- c(2, 1, 4)
  param <- c(N, odds)
  d <- length(odds)

  B <- 4L
  agg_fix <- createMultinomialFamilyCGF(iidReps = B, op_name = "MultinomialCGF_fixedB")

  t_block <- c(0.05, -0.03, 0.01)
  tvec_ok <- rep(t_block, B)
  tvec_bad <- rep(t_block, B - 1L)

  expect_equal(agg_fix$K(tvec_ok, param), B * base$K(t_block, param))
  expect_error(agg_fix$K(tvec_bad, param))  # wrong length should fail
})




test_that("MultinomialCGF aggregates blockwise correctly using MultinomialCGF object", {
  cgf <- MultinomialCGF

  N <- 7
  odds <- c(2, 1, 4)
  param <- c(N, odds)

  d <- length(odds)
  B <- 4L

  t_block <- c(0.05, -0.03, 0.01)
  tvec <- rep(t_block, B)

  # K: sum over blocks
  ref_K <- sum(vapply(seq_len(B), function(i) cgf$K(t_block, param), numeric(1)))
  got_K <- cgf$K(tvec, param)
  expect_equal(as.numeric(got_K), as.numeric(ref_K), tolerance = 1e-12)

  # K1: concatenation
  ref_K1_block <- cgf$K1(t_block, param)
  ref_K1 <- rep(as.numeric(ref_K1_block), B)
  got_K1 <- as.numeric(cgf$K1(tvec, param))
  expect_equal(got_K1, ref_K1, tolerance = 1e-12)

  # K2: block diagonal
  K2_block <- as.matrix(cgf$K2(t_block, param))
  K2_ref <- matrix(0, nrow = d * B, ncol = d * B)
  for (b in seq_len(B)) {
    idx <- ((b - 1L) * d + 1L):(b * d)
    K2_ref[idx, idx] <- K2_block
  }
  got_K2 <- as.matrix(cgf$K2(tvec, param))
  expect_equal(got_K2, K2_ref, tolerance = 1e-12)
})



test_that("MultinomialModelCGF adapts theta and respects prob_vec length", {
  # pick d = 2 (short) or d = 5 (long) to stress dimension handling
  d <- 5L

  # n fixed, prob_vec depends on theta but has fixed length d
  n_fun <- adaptor(fixed_param = 10)
  odds_fun <- function(theta) exp(theta[seq_len(d)])  # positive odds, length d

  cgf <- MultinomialModelCGF(n = n_fun, prob_vec = odds_fun, iidReps = 3)

  theta <- seq(-0.2, 0.2, length.out = d)  # length d
  B <- 3L

  t_block <- seq(-0.1, 0.1, length.out = d)
  tvec <- rep(t_block, B)

  # manual closed form with p = odds / sum(odds)
  odds <- odds_fun(theta)
  p <- odds / sum(odds)
  denom <- sum(p * exp(t_block))
  K_block <- 10 * log(denom)
  K_ref <- B * K_block

  expect_equal(as.numeric(cgf$K(tvec, theta)), as.numeric(K_ref), tolerance = 1e-12)

  # scaling invariance: multiplying odds by a constant must not change K
  odds_fun2 <- function(theta) 7 * exp(theta[seq_len(d)])
  cgf2 <- MultinomialModelCGF(n = n_fun, prob_vec = odds_fun2, iidReps = 3)
  expect_equal(as.numeric(cgf2$K(tvec, theta)), as.numeric(cgf$K(tvec, theta)), tolerance = 1e-12)
})

