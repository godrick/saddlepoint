
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

  harmless_underflow_factor <- function(tvec, p, A) {
    c(
      saddlepoint:::.K2_dense_term(A %*% t(A)),
      saddlepoint:::.K2_factor_term(
        A %*% matrix(1e-200, 1, 1),
        1e-300 * (1 + p[1]^2)
      )
    )
  }
  harmless_underflow_objective <- function(p) {
    saddlepoint:::.K2_factor_logdet(
      harmless_underflow_factor, 0, p
    ) + saddlepoint:::.K2_factor_solve(
      harmless_underflow_factor, 0, p, 2
    )
  }
  harmless_underflow_tape <- RTMB::MakeTape(
    harmless_underflow_objective, 0
  )
  expect_equal(
    c(
      harmless_underflow_tape(0),
      harmless_underflow_tape$jacobian(0),
      harmless_underflow_tape$jacfun()$jacobian(0)
    ),
    c(2, 0, 0),
    tolerance = 0
  )
})

test_that("mapped func_T rejects invalid AD evaluations and recovers", {
  child <- createCGF(
    K = function(tvec, p) {
      0.5 * (p[1] * tvec[1]^2 + tvec[2]^2) +
        (p[1] * tvec[1]^4 + tvec[2]^4) / 24
    },
    K1 = function(tvec, p) {
      c(p[1] * (tvec[1] + tvec[1]^3 / 6),
        tvec[2] + tvec[2]^3 / 6)
    },
    K2 = function(tvec, p) {
      diag(c(p[1] * (1 + tvec[1]^2 / 2),
             1 + tvec[2]^2 / 2), 2L)
    },
    K3operator = function(tvec, p, a, b, c) {
      p[1] * tvec[1] * a[1] * b[1] * c[1] +
        tvec[2] * a[2] * b[2] * c[2]
    },
    K4operator = function(tvec, p, a, b, c, d) {
      p[1] * a[1] * b[1] * c[1] * d[1] +
        a[2] * b[2] * c[2] * d[2]
    },
    K2_factor = function(tvec, p, A) {
      saddlepoint:::.K2_factor_term(
        A,
        c(p[1] * (1 + tvec[1]^2 / 2), 1 + tvec[2]^2 / 2)
      )
    }
  )
  mapped <- linearlyMappedCGF(
    child,
    matrix(c(1, 1), nrow = 1L),
    iidReps = 1L
  )
  tape <- RTMB::MakeTape(
    function(p) mapped$.private_api$func_T(0, p),
    1
  )

  expect_equal(
    c(tape(1), tape$jacobian(1), tape$jacfun()$jacobian(1)),
    c(1 / 16, -1 / 32, 1 / 32),
    tolerance = 1e-12
  )
  expect_true(all(is.nan(c(
    tape(-2), tape$jacobian(-2), tape$jacfun()$jacobian(-2)
  ))))
  expect_equal(
    c(tape(3), tape$jacobian(3), tape$jacfun()$jacobian(3)),
    c(1 / 32, -1 / 128, 1 / 256),
    tolerance = 1e-12
  )
})

test_that("mapped func_T accepts solve roundoff skew and rejects asymmetry", {
  n <- 5L
  set.seed(313)
  directions <- qr.Q(qr(matrix(rnorm(n * n), nrow = n)))
  eigenvalues <- exp(seq(0, -log(3000), length.out = n))
  root <- sweep(directions, 2L, sqrt(eigenvalues), "*")
  map <- cbind(root, 0)

  child <- PoissonModelCGF(
    lambda = function(p) rep(exp(p[1]), n + 1L),
    iidReps = 1L
  )
  mapped <- linearlyMappedCGF(child, map, iidReps = 1L)
  tvec <- numeric(n)
  Q <- mapped$K2_solve(tvec, 0, diag(n))

  expect_gt(
    min(eigen((Q + t(Q)) / 2, symmetric = TRUE, only.values = TRUE)$values),
    0
  )
  expect_silent(chol(Q))
  expect_equal(mapped$.private_api$func_T(tvec, 0), -5 / 12,
               tolerance = 1e-10)

  tape <- RTMB::MakeTape(
    function(p) mapped$.private_api$func_T(tvec, p),
    0
  )
  expect_equal(
    c(tape(0), tape$jacobian(0), tape$jacfun()$jacobian(0)),
    c(-5 / 12, 5 / 12, -5 / 12),
    tolerance = 1e-10
  )

  roundoff_Q <- diag(n)
  roundoff_Q[1L, 2L] <- 0.2 + 1e-13
  roundoff_Q[2L, 1L] <- 0.2 - 1e-13
  symmetric_Q <- (roundoff_Q + t(roundoff_Q)) / 2
  expect_equal(
    mapped$K4operatorAABB(tvec, 0, roundoff_Q),
    mapped$K4operatorAABB(tvec, 0, symmetric_Q),
    tolerance = 1e-12
  )

  asymmetric_Q <- diag(n)
  asymmetric_Q[1L, 2L] <- 0.1
  asymmetric_Q[2L, 1L] <- -0.1
  expect_error(
    mapped$K4operatorAABB(tvec, 0, asymmetric_Q),
    "singular|indefinite|numerically invalid"
  )
})
test_that("solver tvec atomic uses the public K2_solve method", {
  calls <- new.env(parent = emptyenv())
  calls$dense <- 0L
  calls$solve <- 0L
  cgf <- createCGF(
    K = function(tvec, p) 0.5 * exp(p[1]) * sum(tvec^2),
    K1 = function(tvec, p) exp(p[1]) * tvec,
    K2 = function(tvec, p) {
      calls$dense <- calls$dense + 1L
      stop("dense K2 must not be called")
    },
    K3operator = function(tvec, p, a, b, c) 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) 0 * p[1],
    K2_solve = function(tvec, p, rhs) {
      calls$solve <- calls$solve + 1L
      rhs / exp(p[1])
    }
  )
  solver <- function(theta, y, cgf, starting.tvec) y / exp(theta[1])
  t_atomic <- saddlepoint:::make_solver_tvec_atomic(
    cgf, y = 2, solver_fun = solver, theta_init = 0
  )
  tape <- RTMB::MakeTape(function(p) t_atomic(p)[1L], 0)

  expect_equal(
    c(tape(0), tape$jacobian(0), tape$jacfun()$jacobian(0)),
    c(2, -2, 2),
    tolerance = 1e-12
  )
  expect_identical(calls$dense, 0L)
  expect_gt(calls$solve, 0L)
})

test_that("factor terminals validate the completed covariance and scale safely", {
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
    "singular|indefinite|numerically invalid"
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

  multi_rhs <- cbind(rhs, rep(0, 3))
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      scaled_factor,
      numeric(3),
      numeric(),
      multi_rhs
    ),
    cbind(solution, rep(0, 3)),
    tolerance = 1e-10
  )

  # A structurally exact zero solution coordinate remains valid.
  diagonal_root <- diag(c(0.5, 2))
  diagonal_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% diagonal_root, c(1, 1))
  }
  diagonal_rhs <- as.vector(tcrossprod(diagonal_root) %*% c(0, 3))
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      diagonal_factor,
      numeric(2),
      numeric(),
      diagonal_rhs
    ),
    c(0, 3),
    tolerance = 1e-14
  )
})

test_that("dense factors remain finite near the largest double scale", {
  large_scale <- 0.75 * .Machine$double.xmax
  Q <- large_scale * matrix(c(1, 0.25, 0.25, 1), nrow = 2)
  factor <- saddlepoint:::.K2_dense_spd_factor(Q, normalize = FALSE)

  expect_true(all(is.finite(factor$B)))
  expect_equal(
    tcrossprod(factor$B / sqrt(large_scale)),
    Q / large_scale,
    tolerance = 1e-14
  )
})

test_that("splitting factor columns and adding zero columns preserves results", {
  K2 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  root <- t(chol(K2))
  split_count <- 100L
  B <- root[, rep(seq_len(2), each = split_count), drop = FALSE]
  d <- rep(1 / split_count, 2L * split_count)
  B <- cbind(B, matrix(0, nrow = 2, ncol = 3))
  d <- c(d, rep(1, 3))
  factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% B, d)
  }
  rhs <- c(0.2, -0.3)

  expect_equal(
    saddlepoint:::.K2_factor_logdet(factor, numeric(2), numeric()),
    as.numeric(determinant(K2, logarithm = TRUE)$modulus),
    tolerance = 1e-12
  )
  expect_equal(
    saddlepoint:::.K2_factor_solve(factor, numeric(2), numeric(), rhs),
    solve(K2, rhs),
    tolerance = 1e-12
  )
})

test_that("direct-factor QR recovers rank hidden by a rounded Gram matrix", {
  B <- matrix(c(1, 1, 0, 1e-8), nrow = 2)
  row_scale <- sqrt(rowSums(B^2))
  rounded_gram <- tcrossprod(B / row_scale)
  factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% B, c(1, 1))
  }

  expect_error(chol(rounded_gram))
  expect_equal(
    saddlepoint:::.K2_factor_logdet(factor, numeric(2), numeric()),
    2 * log(1e-8),
    tolerance = 1e-12
  )
})

test_that("factor solve applies a componentwise backward-residual sanity check", {
  R <- chol(matrix(c(2, 0.3, 0.3, 1), nrow = 2))
  z <- c(0.5, -1)
  w <- solve(R, solve(t(R), z))

  expect_equal(
    saddlepoint:::.K2_residual_check(w, R, w, z),
    w,
    tolerance = 0
  )

  damaged_w <- w
  damaged_w[1] <- damaged_w[1] + 1e-4
  expect_true(all(is.nan(
    saddlepoint:::.K2_residual_check(damaged_w, R, damaged_w, z)
  )))
})

test_that("rank-completing compositions retain solve and logdet derivatives", {
  count <- PoissonModelCGF(lambda = function(x) x[1], iidReps = 1L)
  summand <- MultinomialModelCGF(
    n = adaptor(fixed_param = 1),
    prob_vec = function(x) x[2:4],
    iidReps = 1L
  )
  rss <- randomlyStoppedSumCGF(
    count, summand, block_size = 3L, iidReps = 1L
  )
  wrapped_summand <- sumOfIndependentCGF(list(summand), iidReps = 1L)
  wrapped_rss <- randomlyStoppedSumCGF(
    count, wrapped_summand, block_size = 3L, iidReps = 1L
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
  expect_equal(wrapped_rss$K2(tvec, theta), rss$K2(tvec, theta))
  expect_equal(wrapped_rss$logdetK2(tvec, theta), rss$logdetK2(tvec, theta))
  expect_equal(
    wrapped_rss$K2_solve(tvec, theta, rhs),
    rss$K2_solve(tvec, theta, rhs),
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

test_that("mapped numeric factor solve remains accurate above the former cutoff", {
  A <- matrix(c(
    -1.0295831616329141, 0.7537973286395325, 0.3180092254628808,
    -1.6734795183709306, 0.3182798250721517, 1.3662313469694245,
    1.6443574037225106, 0.6203839204135480, 1.1389923430543860,
    0.0183791078191809, 2.6853125516886749, 1.2013991211883255,
    -0.8181043471010678, -0.4705873686301221, -2.1426447924667360,
    0.8165912518369038, 0.4653284105581280, 1.3121124773926873,
    -0.2387924084851842, -1.4680435895916326, 0.8322544119216485,
    -0.1533399633773679, -1.1917775586639732, 1.5008761099640489
  ), nrow = 4L)
  parameter <- c(
    9, 123.6878902346769, 0.0807027826860322,
    75.12987031131783, 0.0069805450341088,
    0.1110180815031406, 0.0305517629605729
  )
  tvec <- c(
    0.0202693046070635, 0.0636407276615500,
    -0.0470632877200842, -0.0163980430923402
  )
  rhs <- matrix(c(
    -0.895832957150837, -1.75552178844533,
    0.587442549689941, 0.944209293348568,
    -0.658915585548328, 1.92977712373801,
    -1.86136135212416, -2.28421537253286
  ), nrow = 4L)
  mapped <- linearlyMappedCGF(MultinomialCGF, A, iidReps = 1L)
  K2 <- as.matrix(mapped$K2(tvec, parameter))
  condition <- kappa(K2, exact = TRUE)
  solution <- mapped$K2_solve(tvec, parameter, rhs)

  expect_gt(condition, 9000)
  expect_lt(condition, 9500)
  expect_equal(solution, solve(K2, rhs), tolerance = 1e-10)
  expect_lt(
    max(abs(K2 %*% solution - rhs)) /
      (max(abs(K2)) * max(abs(solution)) + max(abs(rhs))),
    1e-13
  )
})

test_that("numeric factor solve remains available around condition 67825", {
  condition_number <- 67825
  rho <- (condition_number - 1) / (condition_number + 1)
  K2 <- matrix(c(1, rho, rho, 1), nrow = 2)
  root <- t(chol(K2))
  target <- c(1, -0.7)
  rhs <- as.vector(K2 %*% target)
  factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% root, c(1, 1))
  }

  solution <- saddlepoint:::.K2_factor_solve(
    factor, numeric(2), numeric(), rhs
  )
  denominator <- abs(rhs) + abs(K2) %*% abs(solution)
  backward_error <- max(abs(K2 %*% solution - rhs) / denominator)

  expect_true(all(is.finite(solution)))
  expect_equal(as.vector(solution), target, tolerance = 1e-11)
  expect_lt(backward_error, 1e-13)
})

test_that("parameter-dependent taped solves remain available above condition 8192", {
  factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A,
      c(parameter_vector[1], 1)
    )
  }
  rhs <- c(1, 0.25)
  objective <- function(parameter_vector) {
    sum(saddlepoint:::.K2_factor_solve(
      factor, numeric(2), parameter_vector, rhs
    ))
  }
  tape <- RTMB::MakeTape(objective, 0.1)
  derivative_tape <- tape$jacfun()

  for (condition_number in c(8334, 20001, 67825)) {
    p <- 1 / condition_number
    observed <- c(
      tape(p),
      tape$jacobian(p),
      derivative_tape$jacobian(p)
    )
    expected <- c(
      1 / p + 0.25,
      -1 / p^2,
      2 / p^3
    )

    expect_true(all(is.finite(observed)))
    expect_equal(observed / expected, rep(1, 3), tolerance = 1e-12)
  }
})


test_that("mapped multinomial logdet remains accurate across ordinary conditions", {
  selector <- rbind(c(1, 0, 0), c(0, 1, 0))
  mapped <- linearlyMappedCGF(
    MultinomialCGF, selector, iidReps = 1L
  )
  multinomial_parameters <- function(x) {
    parameters <- x[rep(1L, 4L)] * 0
    parameters[1:3] <- 1
    parameters[4] <- exp(x[1])
    parameters
  }

  logdet_objective <- function(x) {
    mapped$logdetK2(numeric(2), multinomial_parameters(x))
  }
  logdet_tape <- RTMB::MakeTape(logdet_objective, 0)
  logdet_derivative_tape <- logdet_tape$jacfun()
  evaluate_logdet <- function(x) {
    c(
      value = logdet_tape(x),
      gradient = logdet_tape$jacobian(x),
      hessian = logdet_derivative_tape$jacobian(x)
    )
  }

  for (condition_number in c(201, 801, 2001, 8334, 20001, 200001)) {
    q <- 2 / (condition_number - 1)
    x <- log(q)
    expected_logdet <- c(
      value = x - 3 * log(2 + q),
      gradient = 1 - 3 * q / (2 + q),
      hessian = -6 * q / (2 + q)^2
    )
    expect_equal(
      evaluate_logdet(x), expected_logdet,
      tolerance = 1e-9
    )
  }
})

test_that("logdet derivatives and invalid evaluations recover", {
  singular_root <- cbind(c(1, 1), c(1, -1))
  singular_factor <- function(tvec, p, A) {
    saddlepoint:::.K2_factor_term(
      A %*% singular_root, c(1 + 0 * p[1], p[1])
    )
  }
  singular_tape <- RTMB::MakeTape(function(p) {
    saddlepoint:::.K2_factor_logdet(
      singular_factor, numeric(2), p
    )
  }, 0.5)
  singular_derivative_tape <- singular_tape$jacfun()
  evaluate_singular <- function(p) {
    c(
      value = singular_tape(p),
      gradient = singular_tape$jacobian(p),
      hessian = singular_derivative_tape$jacobian(p)
    )
  }
  expect_equal(
    evaluate_singular(0.5),
    c(value = log(2), gradient = 2, hessian = -4),
    tolerance = 1e-10
  )
  expect_equal(
    evaluate_singular(1e-6),
    c(value = log(4e-6), gradient = 1e6, hessian = -1e12),
    tolerance = 1e-8
  )
  expect_true(all(is.nan(evaluate_singular(0))))
  expect_equal(
    evaluate_singular(0.25),
    c(value = 0, gradient = 4, hessian = -16),
    tolerance = 1e-10
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

  # Factoring an SPD Q for a mapped higher-order contraction is not a terminal
  # inverse operation.  A valid, near-singular Q must therefore remain usable.
  poisson_child <- PoissonModelCGF(
    lambda = adaptor(fixed_param = c(1, 2, 3)), iidReps = 1L
  )
  reducing_map <- matrix(c(1, 0, 1, 0, 1, -1), 2, byrow = TRUE)
  mapped_poisson <- linearlyMappedCGF(
    poisson_child, reducing_map, iidReps = 1L
  )
  near_singular_Q <- matrix(c(
    1, 1 - 1e-8,
    1 - 1e-8, 1
  ), 2)
  pulled_Q <- t(reducing_map) %*% near_singular_Q %*% reducing_map
  expect_equal(
    mapped_poisson$K4operatorAABB(numeric(2), 0, near_singular_Q),
    sum(c(1, 2, 3) * diag(pulled_Q)^2),
    tolerance = 1e-12
  )
})

test_that("RSS preserves authoritative structured covariance routes", {
  d <- 6L
  diagonal <- seq(1, 2, length.out = d)
  calls <- new.env(parent = emptyenv())
  calls$solve <- 0L
  calls$logdet <- 0L
  child <- createCGF(
    K = function(tvec, p) sum(diagonal * (exp(tvec) - 1)) + 0 * p[1],
    K1 = function(tvec, p) diagonal * exp(tvec) + 0 * p[1],
    K2 = function(tvec, p) diag(diagonal * exp(tvec), d) + 0 * p[1],
    K3operator = function(tvec, p, a, b, c) {
      sum(diagonal * exp(tvec) * a * b * c) + 0 * p[1]
    },
    K4operator = function(tvec, p, a, b, c, z) {
      sum(diagonal * exp(tvec) * a * b * c * z) + 0 * p[1]
    },
    K2_solve = function(tvec, p, rhs) {
      calls$solve <- calls$solve + 1L
      rhs / (diagonal * exp(tvec)) + 0 * p[1]
    },
    logdetK2 = function(tvec, p) {
      calls$logdet <- calls$logdet + 1L
      sum(log(diagonal) + tvec) + 0 * p[1]
    },
    K2_factor = function(tvec, p, A) {
      saddlepoint:::.K2_factor_term(
        t(t(A) * sqrt(diagonal * exp(tvec))), rep(1, d)
      )
    },
    K2_factor_terminal = function() TRUE
  )
  wrapped_children <- list(
    direct = child,
    adapted = adaptCGF(child, function(p) p),
    shifted = shiftedCGF(child, rep(0, d)),
    tilted = ExponentialTiltCGF(child, rep(0, d), iidReps = 1L),
    summed_iid = sumOfiidCGF(child, n = 1),
    replicated = iidReplicatesCGF(child, block_size = d, iidReps = 1L),
    concatenated = concatenationCGF(
      list(child), component_dims = d, iidReps = 1L
    ),
    independent = sumOfIndependentCGF(list(child), iidReps = 1L)
  )
  count <- PoissonModelCGF(
    lambda = adaptor(fixed_param = 1.2), iidReps = 1L
  )
  tvec <- numeric(d)
  rhs <- seq_len(d) / d

  for (wrapped_child in wrapped_children) {
    calls$solve <- 0L
    calls$logdet <- 0L
    rss <- randomlyStoppedSumCGF(
      count, wrapped_child, block_size = d, iidReps = 1L
    )
    K2 <- as.matrix(rss$K2(tvec, 0))
    expect_equal(rss$K2_solve(tvec, 0, rhs), solve(K2, rhs), tolerance = 1e-11)
    expect_equal(
      rss$logdetK2(tvec, 0),
      as.numeric(determinant(K2, logarithm = TRUE)$modulus),
      tolerance = 1e-11
    )
    expect_identical(c(calls$solve, calls$logdet), c(3L, 1L))
  }

  # A one-sided public override remains authoritative when called directly,
  # but does not declare its missing counterpart safe for an enclosing RSS update.
  # The enclosing model must therefore use the completed final factor.
  make_one_sided_child <- function(method_name) {
    one_sided_calls <- new.env(parent = emptyenv())
    one_sided_calls$count <- 0L
    common <- list(
      K = function(tvec, p) sum(diagonal * (exp(tvec) - 1)) + 0 * p[1],
      K1 = function(tvec, p) diagonal * exp(tvec) + 0 * p[1],
      K2 = function(tvec, p) diag(diagonal * exp(tvec), d) + 0 * p[1],
      K3operator = function(tvec, p, a, b, c) 0 * p[1],
      K4operator = function(tvec, p, a, b, c, z) 0 * p[1],
      K2_factor = function(tvec, p, A) {
        saddlepoint:::.K2_factor_term(
          t(t(A) * sqrt(diagonal * exp(tvec))), rep(1, d)
        )
      },
      K2_factor_terminal = function() TRUE
    )
    override <- if (method_name == "K2_solve") {
      list(K2_solve = function(tvec, p, rhs) {
        one_sided_calls$count <- one_sided_calls$count + 1L
        rhs / (diagonal * exp(tvec)) + 0 * p[1]
      })
    } else {
      list(logdetK2 = function(tvec, p) {
        one_sided_calls$count <- one_sided_calls$count + 1L
        sum(log(diagonal) + tvec) + 0 * p[1]
      })
    }
    list(
      cgf = do.call(createCGF, c(common, override)),
      calls = one_sided_calls
    )
  }

  for (method_name in c("K2_solve", "logdetK2")) {
    one_sided <- make_one_sided_child(method_name)
    expect_false(saddlepoint:::.K2_structured_pair_is_safe(one_sided$cgf))
    if (method_name == "K2_solve") {
      expect_equal(
        one_sided$cgf$K2_solve(tvec, 0, rhs),
        rhs / diagonal,
        tolerance = 1e-12
      )
    } else {
      expect_equal(
        one_sided$cgf$logdetK2(tvec, 0),
        sum(log(diagonal)),
        tolerance = 1e-12
      )
    }
    expect_identical(one_sided$calls$count, 1L)
    one_sided$calls$count <- 0L

    rss <- randomlyStoppedSumCGF(
      count, one_sided$cgf, block_size = d, iidReps = 1L
    )
    rss_K2 <- as.matrix(rss$K2(tvec, 0))
    expect_equal(rss$K2_solve(tvec, 0, rhs), solve(rss_K2, rhs),
                 tolerance = 1e-11)
    expect_equal(
      rss$logdetK2(tvec, 0),
      as.numeric(determinant(rss_K2, logarithm = TRUE)$modulus),
      tolerance = 1e-11
    )
    expect_identical(one_sided$calls$count, 0L)
  }

  # The same rule applies when a singleton structural wrapper inherits a safe
  # pair and replaces only one member through `...`: the direct override keeps
  # precedence, but the mixed pair must not be advertised to an enclosing RSS.
  for (method_name in c("K2_solve", "logdetK2")) {
    mixed_calls <- new.env(parent = emptyenv())
    mixed_calls$count <- 0L
    override <- if (method_name == "K2_solve") {
      list(K2_solve = function(tvec, p, rhs) {
        mixed_calls$count <- mixed_calls$count + 1L
        999 * rhs + 0 * p[1]
      })
    } else {
      list(logdetK2 = function(tvec, p) {
        mixed_calls$count <- mixed_calls$count + 1L
        999 + 0 * p[1]
      })
    }
    mixed <- do.call(
      sumOfIndependentCGF,
      c(list(cgf_list = list(child), iidReps = 1L), override)
    )
    expect_false(saddlepoint:::.K2_structured_pair_is_safe(mixed))
    if (method_name == "K2_solve") {
      expect_equal(mixed$K2_solve(tvec, 0, rhs), 999 * rhs)
    } else {
      expect_equal(mixed$logdetK2(tvec, 0), 999)
    }
    expect_identical(mixed_calls$count, 1L)
    mixed_calls$count <- 0L

    rss <- randomlyStoppedSumCGF(
      count, mixed, block_size = d, iidReps = 1L
    )
    rss_K2 <- as.matrix(rss$K2(tvec, 0))
    expect_equal(
      rss$K2_solve(tvec, 0, rhs), solve(rss_K2, rhs),
      tolerance = 1e-11
    )
    expect_equal(
      rss$logdetK2(tvec, 0),
      as.numeric(determinant(rss_K2, logarithm = TRUE)$modulus),
      tolerance = 1e-11
    )
    expect_identical(mixed_calls$count, 0L)
  }

  # Replacing both members supplies a complete direct pair and remains the
  # authoritative structured route.
  replacement_calls <- new.env(parent = emptyenv())
  replacement_calls$solve <- 0L
  replacement_calls$logdet <- 0L
  replaced_pair <- sumOfIndependentCGF(
    list(child),
    iidReps = 1L,
    K2_solve = function(tvec, p, rhs) {
      replacement_calls$solve <- replacement_calls$solve + 1L
      rhs / (diagonal * exp(tvec)) + 0 * p[1]
    },
    logdetK2 = function(tvec, p) {
      replacement_calls$logdet <- replacement_calls$logdet + 1L
      sum(log(diagonal) + tvec) + 0 * p[1]
    }
  )
  expect_true(saddlepoint:::.K2_structured_pair_is_safe(replaced_pair))
  rss <- randomlyStoppedSumCGF(
    count, replaced_pair, block_size = d, iidReps = 1L
  )
  rss_K2 <- as.matrix(rss$K2(tvec, 0))
  expect_equal(rss$K2_solve(tvec, 0, rhs), solve(rss_K2, rhs),
               tolerance = 1e-11)
  expect_equal(
    rss$logdetK2(tvec, 0),
    as.numeric(determinant(rss_K2, logarithm = TRUE)$modulus),
    tolerance = 1e-11
  )
  expect_identical(
    c(replacement_calls$solve, replacement_calls$logdet), c(3L, 1L)
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
