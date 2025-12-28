
# reference for independent Poisson vector
poisson_vec_K <- function(t, lambda_vec) sum(lambda_vec * (exp(t) - 1))
poisson_vec_K1 <- function(t, lambda_vec) lambda_vec * exp(t)
poisson_vec_K2 <- function(t, lambda_vec) diag(lambda_vec * exp(t), nrow = length(lambda_vec))
poisson_vec_K2op <- function(t, lambda_vec, x, y) sum(lambda_vec * exp(t) * x * y)
poisson_vec_K3op <- function(t, lambda_vec, u1, u2, u3) sum(lambda_vec * exp(t) * u1 * u2 * u3)
poisson_vec_K4op <- function(t, lambda_vec, u1, u2, u3, u4) sum(lambda_vec * exp(t) * u1 * u2 * u3 * u4)

test_that("Requires at least one of block_size or iidReps (disambiguation)", {
  # Minimal scalar building blocks
  count <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summ  <- BinomialModelCGF(n = adaptor(fixed_param = 1), p = adaptor(indices = 2), iidReps = 1)

  expect_error(randomlyStoppedSumCGF(count, summ), "supply at least one")
  expect_error(randomlyStoppedSumCGF(count, summ, iidReps = "any", block_size = NULL), "requires a non-NULL")
})

test_that("Poisson count + Bernoulli summand collapses to Poisson(lambda*p) (scalar)", {
  count <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summ  <- BinomialModelCGF(n = adaptor(fixed_param = 1), p = adaptor(indices = 2), iidReps = 1)

  # Disambiguate via block_size only => iidReps becomes "any"
  rss <- randomlyStoppedSumCGF(count, summ, block_size = 1)

  theta <- c(lambda = 2.0, p = 0.3)
  t0 <- 0.12

  lam_eff <- theta[["lambda"]] * theta[["p"]]

  # Closed form: Y ~ Poisson(lam_eff)
  K_ref  <- lam_eff * (exp(t0) - 1)
  K1_ref <- lam_eff * exp(t0)
  K2_ref <- lam_eff * exp(t0)

  expect_equal(as.numeric(rss$K(t0, theta)), as.numeric(K_ref), tolerance = 1e-12)
  expect_equal(as.numeric(rss$K1(t0, theta)), as.numeric(K1_ref), tolerance = 1e-10)
  expect_equal(as.numeric(rss$K2(t0, theta)), as.numeric(K2_ref), tolerance = 1e-10)

  # Higher-order operators (scalar contractions)
  set.seed(1)
  w1 <- rnorm(1); w2 <- rnorm(1); w3 <- rnorm(1); w4 <- rnorm(1)

  K3_ref <- lam_eff * exp(t0) * w1 * w2 * w3
  K4_ref <- lam_eff * exp(t0) * w1 * w2 * w3 * w4

  expect_equal(as.numeric(rss$K3operator(t0, theta, w1, w2, w3)),
               as.numeric(K3_ref), tolerance = 1e-8)
  expect_equal(as.numeric(rss$K4operator(t0, theta, w1, w2, w3, w4)),
               as.numeric(K4_ref), tolerance = 1e-7)

  # K2_solve + logdetK2 (scalar)
  K2_val <- as.numeric(rss$K2(t0, theta))
  rhs <- 3.7
  expect_equal(as.numeric(rss$K2_solve(t0, theta, rhs)), rhs / K2_val, tolerance = 1e-10)
  expect_equal(as.numeric(rss$logdetK2(t0, theta)), log(K2_val), tolerance = 1e-10)
})

test_that("Replication: block_size only behaves like iidReps='any' for Y (scalar)", {
  count <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summ  <- BinomialModelCGF(n = adaptor(fixed_param = 2), p = adaptor(indices = 2), iidReps = 1)

  rss_any <- randomlyStoppedSumCGF(count, summ, block_size = 1) # iidReps inferred from tvec length
  rss_B   <- randomlyStoppedSumCGF(count, summ, iidReps = 5)    # block_size inferred from tvec length

  theta <- c(lambda = 1.3, p = 0.25)
  B <- 5L
  t_block <- 0.05
  tvec <- rep(t_block, B)

  # Both should interpret tvec as B independent scalar observations
  K_any <- as.numeric(rss_any$K(tvec, theta))
  K_B   <- as.numeric(rss_B$K(tvec, theta))
  expect_equal(K_any, K_B, tolerance = 1e-12)

  # And both should equal sum over blocks of the 1-block evaluation
  ref <- B * as.numeric(rss_any$K(t_block, theta))
  expect_equal(K_any, ref, tolerance = 1e-12)

  # K1 concatenates
  ref_K1 <- rep(as.numeric(rss_any$K1(t_block, theta)), B)
  expect_equal(as.numeric(rss_any$K1(tvec, theta)), ref_K1, tolerance = 1e-12)
  expect_equal(as.numeric(rss_B$K1(tvec, theta)),   ref_K1, tolerance = 1e-12)

  # K2 is block diagonal (for scalar this is just diagonal)
  diag_ref <- rep(as.numeric(rss_any$K2(t_block, theta)), B)
  expect_equal(diag(as.matrix(rss_any$K2(tvec, theta))), diag_ref, tolerance = 1e-12)
  expect_equal(diag(as.matrix(rss_B$K2(tvec, theta))),   diag_ref, tolerance = 1e-12)
})

test_that("Poisson count + categorical summand (multinomial n=1) gives independent Poisson vector", {
  d <- 3L

  count <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summ  <- MultinomialModelCGF(
    n        = adaptor(fixed_param = 1),
    prob_vec = adaptor(indices = 2:(1+d)),
    iidReps  = 1
  )

  rss <- randomlyStoppedSumCGF(count, summ, block_size = d)  # iidReps="any" implied

  theta <- c(lambda = 2.4, p1 = 0.2, p2 = 0.3, p3 = 0.5)
  p <- as.numeric(theta[2:(1+d)])
  lambda_vec <- theta[["lambda"]] * p

  t0 <- c(0.10, -0.20, 0.05)

  expect_equal(as.numeric(rss$K(t0, theta)),
               as.numeric(poisson_vec_K(t0, lambda_vec)),
               tolerance = 1e-10)

  expect_equal(as.numeric(rss$K1(t0, theta)),
               as.numeric(poisson_vec_K1(t0, lambda_vec)),
               tolerance = 1e-8)

  expect_equal(as.matrix(rss$K2(t0, theta)),
               as.matrix(poisson_vec_K2(t0, lambda_vec)),
               tolerance = 1e-8)

  set.seed(2)
  u1 <- rnorm(d); u2 <- rnorm(d); u3 <- rnorm(d); u4 <- rnorm(d)
  expect_equal(as.numeric(rss$K2operator(t0, theta, u1, u2)),
               as.numeric(poisson_vec_K2op(t0, lambda_vec, u1, u2)),
               tolerance = 1e-8)
  expect_equal(as.numeric(rss$K3operator(t0, theta, u1, u2, u3)),
               as.numeric(poisson_vec_K3op(t0, lambda_vec, u1, u2, u3)),
               tolerance = 1e-7)
  expect_equal(as.numeric(rss$K4operator(t0, theta, u1, u2, u3, u4)),
               as.numeric(poisson_vec_K4op(t0, lambda_vec, u1, u2, u3, u4)),
               tolerance = 1e-6)

  # Replication over B blocks
  B <- 4L
  tvec <- rep(t0, B)
  refK <- B * poisson_vec_K(t0, lambda_vec)
  expect_equal(as.numeric(rss$K(tvec, theta)), as.numeric(refK), tolerance = 1e-10)

  # K2 block diagonal: diagonal blocks identical
  K2_block <- poisson_vec_K2(t0, lambda_vec)
  K2_ref <- matrix(0, nrow = d*B, ncol = d*B)
  for (b in seq_len(B)) {
    idx <- ((b-1L)*d + 1L):(b*d)
    K2_ref[idx, idx] <- K2_block
  }
  expect_equal(as.matrix(rss$K2(tvec, theta)), K2_ref, tolerance = 1e-8)
})
