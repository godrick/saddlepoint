test_that("RSS fast K2_solve and logdetK2 agree with explicit solve/determinant", {
  # Use a small vector-valued RSS where K2 is PD.
  # Count: Poisson(lambda_N)
  # Summand: independent Poisson vector (lambda_X1, lambda_X2, lambda_X3)
  #
  d <- 3

  count_cgf <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summand_cgf <- PoissonModelCGF(lambda = adaptor(indices = 2:(d + 1)), iidReps = "any")

  rss <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf,
    block_size = d,
    iidReps = 1
  )

  theta <- c(1.5, 0.7, 1.1, 0.9)  # lambda_N, lambda_X1:lambda_X3
  tvec <- c(0.05, -0.02, 0.03)

  K2 <- as.matrix(rss$K2(tvec, theta))
  expect_true(all(is.finite(K2)))

  # Vector rhs
  set.seed(1)
  rhs <- rnorm(d)
  x_fast <- rss$K2_solve(tvec, theta, rhs)
  x_ref  <- solve(K2, rhs)
  expect_equal(as.numeric(x_fast), as.numeric(x_ref), tolerance = 1e-10)

  # Matrix rhs (multiple right-hand sides)
  RHS <- matrix(rnorm(d * 2), nrow = d, ncol = 2)
  X_fast <- rss$K2_solve(tvec, theta, RHS)
  X_ref  <- solve(K2, RHS)
  expect_equal(as.matrix(X_fast), as.matrix(X_ref), tolerance = 1e-10)

  # logdet
  ld_fast <- as.numeric(rss$logdetK2(tvec, theta))
  ld_ref  <- as.numeric(determinant(K2, logarithm = TRUE)$modulus)
  expect_equal(ld_fast, ld_ref, tolerance = 1e-10)
})



test_that("RSS fast AABB / K3K3 operators and func_T agree with base factored implementations (small d)", {
  # The RSS implementation provides closed forms for:
  #   - K4operatorAABB(Q,Q)
  #   - K3K3operatorAABBCC(Q,Q,Q)
  #   - K3K3operatorABCABC(Q,Q,Q)
  #   - func_T
  #
  # Here we compare them against the base-class *factored* contractions, which are slow but generic.

  d <- 3
  count_cgf <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summand_cgf <- PoissonModelCGF(lambda = adaptor(indices = 2:(d + 1)), iidReps = "any")

  rss <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf,
    block_size = d,
    iidReps = 1
  )

  theta <- c(1.2, 0.6, 1.0, 0.8)
  tvec  <- c(0.04, -0.03, 0.02)

  # Q is the inverse Hessian at tvec
  K2 <- as.matrix(rss$K2(tvec, theta))
  Q  <- solve(K2)
  Q  <- 0.5 * (Q + t(Q))  # symmetrize ??? not necessary it shoulb be already symmetric

  # Build the (A,d) factorization used by the base factored methods:
  #   Q = A diag(d) A^T
  # with A columns scaled from chol(Q).
  U <- chol(Q)            # upper triangular, Q = t(U) %*% U
  diagU <- diag(U)
  dvec <- diagU^2
  A <- t(U) %*% diag(1 / diagU)

  # Base factored methods from CGF class
  K4_factored    <- rss$.private_api$K4operatorAABB_factored
  K3K3_AABB_fact <- rss$.private_api$K3K3operatorAABBCC_factored
  K3K3_ABC_fact  <- rss$.private_api$K3K3operatorABCABC_factored

  # Reference  values
  K4_ref    <- K4_factored(tvec, theta, A, dvec, A, dvec)
  AABB_ref  <- K3K3_AABB_fact(tvec, theta, A, dvec, A, dvec, A, dvec)
  ABC_ref   <- K3K3_ABC_fact(tvec, theta, A, dvec, A, dvec, A, dvec)

  # new methods
  K4_fast   <- rss$K4operatorAABB(tvec, theta, Q, Q)
  AABB_fast <- rss$K3K3operatorAABBCC(tvec, theta, Q, Q, Q)
  ABC_fast  <- rss$K3K3operatorABCABC(tvec, theta, Q, Q, Q)

  expect_equal(as.numeric(K4_fast),  as.numeric(K4_ref),   tolerance = 1e-8)
  expect_equal(as.numeric(AABB_fast), as.numeric(AABB_ref), tolerance = 1e-8)
  expect_equal(as.numeric(ABC_fast),  as.numeric(ABC_ref),  tolerance = 1e-8)

  # func_T comparison
  func_T_fast <- rss$.private_api$func_T
  T_fast <- as.numeric(func_T_fast(tvec, theta))
  T_ref  <- as.numeric(K4_ref / 8 - AABB_ref / 8 - ABC_ref / 12)
  expect_equal(T_fast, T_ref, tolerance = 1e-8)
})
