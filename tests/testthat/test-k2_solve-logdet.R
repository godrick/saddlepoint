
test_that("K2_solve and logdetK2 agree with solve and determinant", {
  cgf <- MultinomialCGF

  # Use a mapped multinomial to make K2 invertible
  A <- matrix(c(1, 0, 0,
                0, 1, 0), nrow = 2, byrow = TRUE)
  mapped <- linearlyMappedCGF(cgf = cgf, matrix_A = A, iidReps = "any")

  theta <- c(12, 2, 3, 5)
  tvec <- c(0.10, -0.05)

  K2 <- as.matrix(mapped$K2(tvec, theta))
  rhs <- c(0.2, -0.3)

  sol_ref <- solve(K2, rhs)
  sol_got <- mapped$K2_solve(tvec, theta, rhs)
  expect_equal(as.numeric(sol_got), as.numeric(sol_ref), tolerance = 1e-8)

  ld_ref <- as.numeric(determinant(K2, logarithm = TRUE)$modulus)
  ld_got <- as.numeric(mapped$logdetK2(tvec, theta))
  expect_equal(ld_got, ld_ref, tolerance = 1e-8)

  summand <- SubunitaryMultinomialModelCGF(
    n = adaptor(fixed_param = 10),
    prob_vec = adaptor(fixed_param = c(0.2, 0.3, 0.1)),
    iidReps = 1L
  )
  dynamic_map <- function(x) {
    out <- matrix(0, nrow = 2, ncol = 3) * x[1]
    out[1, 1] <- 1
    out[2, 2] <- 1
    out[1, 3] <- x[1]
    out[2, 3] <- 0.2 * x[1]
    out
  }
  dynamic <- linearlyMappedCGF(summand, dynamic_map, iidReps = 1L)
  rhs <- cbind(c(1, -1), c(0.5, 0.25))
  theta <- 0.4
  K2 <- as.matrix(dynamic$K2(tvec, theta))
  expect_equal(
    dynamic$logdetK2(tvec, theta),
    as.numeric(determinant(K2, logarithm = TRUE)$modulus),
    tolerance = 1e-12
  )
  expect_equal(
    dynamic$K2_solve(tvec, theta, rhs),
    solve(K2, rhs),
    tolerance = 1e-12
  )

  objective <- function(x) {
    dynamic$logdetK2(tvec, x) +
      sum(dynamic$K2_solve(tvec, x, rhs))
  }
  tape <- RTMB::MakeTape(objective, theta)
  expect_equal(
    as.numeric(tape$jacobian(theta)),
    as.numeric(numDeriv::grad(objective, theta)),
    tolerance = 1e-7
  )
  expect_equal(
    as.numeric(tape$jacfun()$jacobian(theta)),
    as.numeric(numDeriv::hessian(objective, theta)),
    tolerance = 1e-6
  )
})

test_that("factor terminals preserve zero-weight derivatives and recover", {
  b <- c(1, 2)
  rhs <- c(0.5, -1)
  make_cgf <- function(weight) {
    createCGF(
      K = function(tvec, parameter_vector) 0 * parameter_vector[1],
      K1 = function(tvec, parameter_vector) 0 * tvec,
      K2 = function(tvec, parameter_vector) {
        diag(2) + weight(parameter_vector) * tcrossprod(b)
      },
      K3operator = function(tvec, parameter_vector, a, b, c) {
        0 * parameter_vector[1]
      },
      K4operator = function(tvec, parameter_vector, a, b, c, d) {
        0 * parameter_vector[1]
      },
      K2_factor = function(tvec, parameter_vector, A) {
        saddlepoint:::.K2_factor_term(
          A %*% cbind(b, diag(2)),
          c(weight(parameter_vector), 1, 1)
        )
      },
      K2_factor_terminal = function() TRUE
    )
  }

  cgf <- make_cgf(function(x) x[1])
  objective <- function(x) {
    cgf$logdetK2(numeric(2), x) +
      sum(cgf$K2_solve(numeric(2), x, rhs))
  }
  tape <- RTMB::MakeTape(objective, 0)
  expect_equal(as.numeric(tape(0)), -0.5, tolerance = 1e-12)
  expect_equal(as.numeric(tape$jacobian(0)), 9.5, tolerance = 1e-10)
  expect_equal(
    as.numeric(tape$jacfun()$jacobian(0)),
    -70,
    tolerance = 1e-8
  )

  expect_true(is.nan(as.numeric(tape(-0.1))))
  expect_true(is.nan(as.numeric(tape$jacobian(-0.1))))
  expect_true(is.nan(as.numeric(tape$jacfun()$jacobian(-0.1))))
  expect_equal(as.numeric(tape(0.2)), log(2) - 0.05, tolerance = 1e-12)
  expect_equal(as.numeric(tape$jacobian(0.2)), 3.625, tolerance = 1e-10)
  expect_equal(
    as.numeric(tape$jacfun()$jacobian(0.2)),
    -11.875,
    tolerance = 1e-8
  )

  underflow_cgf <- make_cgf(function(x) exp(x[1] - 1000))
  underflow_tape <- RTMB::MakeTape(function(x) {
    underflow_cgf$logdetK2(numeric(2), x) +
      sum(underflow_cgf$K2_solve(numeric(2), x, rhs))
  }, 0)
  expect_equal(as.numeric(underflow_tape(0)), -0.5, tolerance = 1e-12)
  expect_equal(as.numeric(underflow_tape$jacobian(0)), 0, tolerance = 0)
  expect_equal(
    as.numeric(underflow_tape$jacfun()$jacobian(0)),
    0,
    tolerance = 0
  )
})

test_that("factor terminals validate completion, scale, and numerical loss", {
  S <- diag(c(-2e-15, rep(0, 29)))
  K2 <- S + diag(30)
  rhs <- seq_len(30) / 30
  factor <- function(tvec, parameter_vector, A) {
    c(
      saddlepoint:::.K2_dense_term(A %*% S %*% t(A)),
      saddlepoint:::.K2_factor_term(A, rep(1, 30))
    )
  }

  expect_equal(min(diag(S)), -2e-15, tolerance = 0)
  expect_gt(min(eigen(K2, symmetric = TRUE, only.values = TRUE)$values), 0.99)
  expect_equal(
    saddlepoint:::.K2_factor_solve(factor, numeric(30), numeric(), rhs),
    solve(K2, rhs),
    tolerance = 2e-12
  )
  expect_equal(
    saddlepoint:::.K2_factor_logdet(factor, numeric(30), numeric()),
    as.numeric(determinant(K2, logarithm = TRUE)$modulus),
    tolerance = 2e-12
  )

  singular_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% matrix(1, 2, 1), 1)
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      singular_factor,
      numeric(2),
      numeric(),
      c(1, 1)
    ),
    "singular|ill-conditioned"
  )

  scales <- c(1e-100, 1, 1e100)
  G <- matrix(c(
    1, 0.2, 0.1,
    0.3, 1, 0.2,
    0.1, 0.4, 1
  ), 3, byrow = TRUE)
  B <- diag(scales) %*% G
  expected <- c(0.5, -1, 2)
  rhs <- scales * (G %*% (t(G) %*% expected))
  scaled_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% B, rep(1, 3))
  }

  solution <- saddlepoint:::.K2_factor_solve(
    scaled_factor,
    numeric(3),
    numeric(),
    rhs
  )
  expect_equal(as.vector(solution * scales), expected, tolerance = 1e-10)
  expect_equal(
    saddlepoint:::.K2_factor_logdet(
      scaled_factor,
      numeric(3),
      numeric()
    ),
    2 * log(abs(det(G))),
    tolerance = 1e-11
  )

  n <- 36L
  R0 <- diag(n)
  R0[cbind(1:35, 2:36)] <- 2
  ill_conditioned_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% t(R0), rep(1, n))
  }
  set.seed(6907)
  x0 <- rnorm(n)
  bad_rhs <- crossprod(R0, R0 %*% x0)
  expect_error(
    saddlepoint:::.K2_factor_solve(
      ill_conditioned_factor,
      numeric(n),
      numeric(),
      bad_rhs
    ),
    "numerical accuracy|numerical loss"
  )
  expect_error(
    saddlepoint:::.K2_factor_logdet(
      ill_conditioned_factor,
      numeric(n),
      numeric()
    ),
    "log-determinant derivatives.*numerical accuracy"
  )

  count <- PoissonModelCGF(lambda = function(x) x[1], iidReps = 1L)
  summand <- MultinomialModelCGF(
    n = adaptor(fixed_param = 1),
    prob_vec = function(x) x[2:4],
    iidReps = 1L
  )
  rss <- randomlyStoppedSumCGF(
    count, summand, block_size = 3L, iidReps = 1L
  )
  theta <- c(2, 0.5, 0.3, 0.2)
  tvec <- rep(0, 3)
  rhs <- c(1, -2, 3)
  expect_equal(
    as.matrix(rss$K2(tvec, theta)),
    diag(c(1, 0.6, 0.4)),
    tolerance = 1e-12
  )
  expect_equal(rss$logdetK2(tvec, theta), log(0.24), tolerance = 1e-12)
  expect_equal(
    rss$K2_solve(tvec, theta, rhs),
    c(1, -10 / 3, 7.5),
    tolerance = 1e-11
  )

  rss_objective <- function(x) {
    rss$logdetK2(tvec, x) + sum(rss$K2_solve(tvec, x, rhs))
  }
  rss_tape <- RTMB::MakeTape(rss_objective, theta)
  expect_equal(
    as.vector(rss_tape$jacobian(theta)),
    numDeriv::grad(rss_objective, theta),
    tolerance = 1e-6
  )
  expect_equal(
    rss_tape$jacfun()$jacobian(theta),
    numDeriv::hessian(rss_objective, theta),
    tolerance = 1e-5
  )
})

test_that("composition order and mixed covariance terms retain one result", {
  X1 <- MultinomialModelCGF(
    n = adaptor(fixed_param = 3),
    prob_vec = adaptor(fixed_param = c(0.2, 0.3, 0.5)),
    iidReps = 1L
  )
  X2 <- MultinomialModelCGF(
    n = adaptor(fixed_param = 4),
    prob_vec = adaptor(fixed_param = c(0.4, 0.1, 0.5)),
    iidReps = 1L
  )
  A <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE)
  map_after_sum <- linearlyMappedCGF(
    sumOfIndependentCGF(list(X1, X2), iidReps = 1L),
    A,
    iidReps = 1L
  )
  sum_after_map <- sumOfIndependentCGF(list(
    linearlyMappedCGF(X1, A, iidReps = 1L),
    linearlyMappedCGF(X2, A, iidReps = 1L)
  ), iidReps = 1L)
  tvec <- c(0.1, -0.2)
  rhs <- c(2, -1)

  expect_equal(map_after_sum$K2(tvec, 1), sum_after_map$K2(tvec, 1))
  expect_equal(
    map_after_sum$logdetK2(tvec, 1),
    sum_after_map$logdetK2(tvec, 1),
    tolerance = 1e-12
  )
  expect_equal(
    map_after_sum$K2_solve(tvec, 1, rhs),
    sum_after_map$K2_solve(tvec, 1, rhs),
    tolerance = 1e-11
  )
  expect_equal(
    map_after_sum$.private_api$func_T(tvec, 1),
    sum_after_map$.private_api$func_T(tvec, 1),
    tolerance = 1e-11
  )

  dense <- createCGF(
    K = function(tvec, p) 0.5 * p[1] * sum(tvec * tvec),
    K1 = function(tvec, p) p[1] * tvec,
    K2 = function(tvec, p) diag(p[1], length(tvec)),
    K3operator = function(tvec, p, a, b, c) 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) 0 * p[1]
  )
  mixed <- sumOfIndependentCGF(list(
    linearlyMappedCGF(X1, A, iidReps = 1L),
    dense
  ), iidReps = 1L)
  theta <- 0.7

  candidate_numerics <- RTMB::MakeTape(function(x) {
    mixed$logdetK2(tvec, x) + sum(mixed$K2_solve(tvec, x, rhs))
  }, theta)
  reference_numerics <- RTMB::MakeTape(function(x) {
    K2 <- mixed$K2(tvec, x)
    determinant(K2, logarithm = TRUE)$modulus + sum(solve(K2, rhs))
  }, theta)
  expect_equal(
    c(
      candidate_numerics(theta),
      candidate_numerics$jacobian(theta),
      candidate_numerics$jacfun()$jacobian(theta)
    ),
    c(
      reference_numerics(theta),
      reference_numerics$jacobian(theta),
      reference_numerics$jacfun()$jacobian(theta)
    ),
    tolerance = 1e-9
  )

  candidate_T <- RTMB::MakeTape(
    function(x) mixed$.private_api$func_T(tvec, x), theta
  )
  reference_T <- RTMB::MakeTape(function(x) {
    Q <- solve(mixed$K2(tvec, x))
    mixed$K4operatorAABB(tvec, x, Q) / 8 -
      mixed$K3K3operatorAABBCC(tvec, x, Q) / 8 -
      mixed$K3K3operatorABCABC(tvec, x, Q) / 12
  }, theta)
  expect_equal(
    c(
      candidate_T(theta),
      candidate_T$jacobian(theta),
      candidate_T$jacfun()$jacobian(theta)
    ),
    c(
      reference_T(theta),
      reference_T$jacobian(theta),
      reference_T$jacfun()$jacobian(theta)
    ),
    tolerance = 1e-9
  )
})

test_that("base func_T uses K2_solve and factored contractions", {
  calls <- new.env(parent = emptyenv())
  calls$solve <- calls$k4 <- calls$k3a <- calls$k3b <- 0L
  cgf <- createCGF(
    K = function(tvec, parameter_vector) 0 * parameter_vector[1],
    K1 = function(tvec, parameter_vector) numeric(length(tvec)) * parameter_vector[1],
    K2 = function(tvec, parameter_vector) diag(length(tvec)) + 0 * parameter_vector[1],
    K3operator = function(tvec, parameter_vector, a, b, c) 0 * parameter_vector[1],
    K4operator = function(tvec, parameter_vector, a, b, c, d) 0 * parameter_vector[1],
    K2_solve = function(tvec, parameter_vector, rhs) {
      calls$solve <- calls$solve + 1L
      rhs + 0 * parameter_vector[1]
    },
    K4operatorAABB_factored = function(tvec, parameter_vector, A, d) {
      calls$k4 <- calls$k4 + 1L
      8 + 0 * parameter_vector[1]
    },
    K3K3operatorAABBCC_factored = function(tvec, parameter_vector, A, d) {
      calls$k3a <- calls$k3a + 1L
      8 + 0 * parameter_vector[1]
    },
    K3K3operatorABCABC_factored = function(tvec, parameter_vector, A, d) {
      calls$k3b <- calls$k3b + 1L
      12 + 0 * parameter_vector[1]
    }
  )

  expect_equal(cgf$.private_api$func_T(c(0, 0), 1), -1)
  expect_identical(
    c(calls$solve, calls$k4, calls$k3a, calls$k3b),
    c(1L, 1L, 1L, 1L)
  )
})
