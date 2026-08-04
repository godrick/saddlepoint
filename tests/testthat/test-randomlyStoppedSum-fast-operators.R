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
  K4_ref    <- K4_factored(tvec, theta, A, dvec)
  AABB_ref  <- K3K3_AABB_fact(tvec, theta, A, dvec)
  ABC_ref   <- K3K3_ABC_fact(tvec, theta, A, dvec)

  # new methods
  K4_fast   <- rss$K4operatorAABB(tvec, theta, Q)
  AABB_fast <- rss$K3K3operatorAABBCC(tvec, theta, Q)
  ABC_fast  <- rss$K3K3operatorABCABC(tvec, theta, Q)

  expect_equal(as.numeric(K4_fast),  as.numeric(K4_ref),   tolerance = 1e-8)
  expect_equal(as.numeric(AABB_fast), as.numeric(AABB_ref), tolerance = 1e-8)
  expect_equal(as.numeric(ABC_fast),  as.numeric(ABC_ref),  tolerance = 1e-8)

  # func_T comparison
  func_T_fast <- rss$.private_api$func_T
  T_fast <- as.numeric(func_T_fast(tvec, theta))
  T_ref  <- as.numeric(K4_ref / 8 - AABB_ref / 8 - ABC_ref / 12)
  expect_equal(T_fast, T_ref, tolerance = 1e-8)
})

test_that("RSS correction dispatch preserves override precedence", {
  count <- PoissonModelCGF(
    lambda = adaptor(fixed_param = 1.2), iidReps = 1L
  )
  summand <- PoissonModelCGF(
    lambda = adaptor(fixed_param = c(0.7, 0.9)), iidReps = 1L
  )
  tvec <- c(0.04, -0.03)
  theta <- 0
  expected <- c(value = -1, gradient = 1 / 8, hessian = 7 / 24)

  make_calls <- function() {
    out <- new.env(parent = emptyenv())
    out$k4 <- out$k3a <- out$k3b <- 0L
    out
  }
  counts <- function(x) c(x$k4, x$k3a, x$k3b)
  values <- list(
    k4 = function(p) 8 + 2 * p[1] + p[1]^2,
    k3a = function(p) 8 - p[1] + p[1]^2 / 2,
    k3b = function(p) 12 + 3 * p[1] - p[1]^2
  )
  dense_methods <- function(calls) list(
    K4operatorAABB = function(tvec, p, Q) {
      calls$k4 <- calls$k4 + 1L
      values$k4(p)
    },
    K3K3operatorAABBCC = function(tvec, p, Q) {
      calls$k3a <- calls$k3a + 1L
      values$k3a(p)
    },
    K3K3operatorABCABC = function(tvec, p, Q) {
      calls$k3b <- calls$k3b + 1L
      values$k3b(p)
    }
  )
  factored_methods <- function(calls) list(
    K4operatorAABB_factored = function(tvec, p, A, d) {
      calls$k4 <- calls$k4 + 1L
      values$k4(p)
    },
    K3K3operatorAABBCC_factored = function(tvec, p, A, d) {
      calls$k3a <- calls$k3a + 1L
      values$k3a(p)
    },
    K3K3operatorABCABC_factored = function(tvec, p, A, d) {
      calls$k3b <- calls$k3b + 1L
      values$k3b(p)
    }
  )
  make_rss <- function(methods = list()) {
    do.call(
      randomlyStoppedSumCGF,
      c(list(count_cgf = count, summand_cgf = summand, iidReps = 1L), methods)
    )
  }
  correction_vgh <- function(cgf) {
    tape <- RTMB::MakeTape(
      function(p) cgf$.private_api$func_T(tvec, p), theta
    )
    c(
      value = tape(theta),
      gradient = tape$jacobian(theta),
      hessian = tape$jacfun()$jacobian(theta)
    )
  }

  for (factory in list(dense_methods, factored_methods)) {
    calls <- make_calls()
    expect_equal(correction_vgh(make_rss(factory(calls))), expected, tolerance = 1e-11)
    expect_true(all(counts(calls) > 0L))
  }

  dense_calls <- make_calls()
  factored_calls <- make_calls()
  both <- make_rss(c(dense_methods(dense_calls), factored_methods(factored_calls)))
  expect_equal(correction_vgh(both), expected, tolerance = 1e-11)
  expect_identical(counts(dense_calls), c(0L, 0L, 0L))
  expect_true(all(counts(factored_calls) > 0L))

  explicit_calls <- 0L
  explicit <- make_rss(list(func_T = function(tvec, p) {
    explicit_calls <<- explicit_calls + 1L
    77 + p[1]^2
  }))
  expect_equal(unname(correction_vgh(explicit)), c(77, 0, 2), tolerance = 1e-12)
  expect_gt(explicit_calls, 0L)
})

test_that("RSS and base K4 bound high-rank factorizations", {
  calls <- new.env(parent = emptyenv())
  calls$k3 <- calls$k4 <- 0L
  count <- PoissonModelCGF(
    lambda = function(p) exp(p[1]), iidReps = 1L
  )
  summand <- createCGF(
    K = function(tvec, p) sum(exp(p[1]) * (exp(tvec) - 1)),
    K1 = function(tvec, p) exp(p[1] + tvec),
    K2 = function(tvec, p) diag(exp(p[1] + tvec), length(tvec)),
    K3operator = function(tvec, p, a, b, c) {
      calls$k3 <- calls$k3 + 1L
      sum(exp(p[1] + tvec) * a * b * c)
    },
    K4operator = function(tvec, p, a, b, c, d) {
      calls$k4 <- calls$k4 + 1L
      sum(exp(p[1] + tvec) * a * b * c * d)
    }
  )
  rss <- randomlyStoppedSumCGF(
    count, summand, block_size = 2L, iidReps = 1L
  )
  set.seed(8128)
  B <- matrix(rnorm(2L * 32L), 2L, 32L) / sqrt(32)
  dvec <- runif(32L, 0.5, 1.5)
  Q <- B %*% (dvec * t(B))
  tvec <- c(0.03, -0.02)
  theta <- 0.1

  reference <- RTMB::MakeTape(function(p) {
    rss$K4operatorAABB(tvec, p, Q)
  }, theta)
  calls$k3 <- calls$k4 <- 0L
  candidate <- RTMB::MakeTape(function(p) {
    rss$.private_api$K4operatorAABB_factored(tvec, p, B, dvec)
  }, theta)
  candidate_calls <- c(k3 = calls$k3, k4 = calls$k4)
  expect_equal(
    c(candidate(theta), candidate$jacobian(theta), candidate$jacfun()$jacobian(theta)),
    c(reference(theta), reference$jacobian(theta), reference$jacfun()$jacobian(theta)),
    tolerance = 1e-11
  )
  expect_lte(unname(candidate_calls[["k4"]]), 8L)
  expect_lte(unname(candidate_calls[["k3"]]), 4L)

  direction <- c(1.2, -0.7)
  rank_deficient_B <- cbind(direction, 2 * direction, -3 * direction)
  rank_deficient_d <- c(0.4, 0.1, 0.05)
  equivalent_d <- sum(rank_deficient_d * c(1, 2, -3)^2)
  wide <- RTMB::MakeTape(function(p) {
    summand$.private_api$K4operatorAABB_factored(
      tvec, p, rank_deficient_B, rank_deficient_d
    )
  }, theta)
  thin <- RTMB::MakeTape(function(p) {
    summand$.private_api$K4operatorAABB_factored(
      tvec, p, matrix(direction, ncol = 1L), equivalent_d
    )
  }, theta)
  expect_equal(
    c(wide(theta), wide$jacobian(theta), wide$jacfun()$jacobian(theta)),
    c(thin(theta), thin$jacobian(theta), thin$jacfun()$jacobian(theta)),
    tolerance = 1e-12
  )
})
