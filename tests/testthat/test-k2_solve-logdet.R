
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

test_that("factor terminals validate only the completed covariance", {
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
})

test_that("factor terminals equilibrate scales and reject numerical loss", {
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
