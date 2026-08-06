
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

  scale_exponent <- 11
  row_scale <- diag(c(10^-scale_exponent, 10^scale_exponent))
  correlation <- matrix(c(1, 0.3, 0.3, 1), 2L)
  unstable_root <- row_scale %*% t(chol(correlation))
  unstable_K2 <- tcrossprod(unstable_root)
  unstable_rhs <- as.vector(unstable_K2 %*% c(1, -0.7))
  unstable_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% unstable_root,
      c(1, 1)
    )
  }

  expect_equal(
    saddlepoint:::.K2_factor_solve(
      unstable_factor,
      numeric(2),
      numeric(),
      c(0, 0)
    ),
    c(0, 0),
    tolerance = 0
  )
  expect_error(
    saddlepoint:::.K2_factor_solve(
      unstable_factor,
      numeric(2),
      numeric(),
      unstable_rhs
    ),
    "numerical accuracy certification"
  )
  expect_error(
    saddlepoint:::.K2_factor_solve(
      unstable_factor,
      numeric(2),
      numeric(),
      cbind(c(0, 0), unstable_rhs)
    ),
    "numerical accuracy certification"
  )

  # A large, accurately solved coordinate must not mask complete loss in a
  # different coordinate, even when the equilibrating row scales are equal.
  correlation_root <- t(chol(correlation))
  correlation_K2 <- tcrossprod(correlation_root)
  no_scale_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% correlation_root,
      c(1, 1)
    )
  }
  no_scale_bad_rhs <- as.vector(correlation_K2 %*% c(0, 1e30))
  expect_error(
    saddlepoint:::.K2_factor_solve(
      no_scale_factor,
      numeric(2),
      numeric(),
      no_scale_bad_rhs
    ),
    "numerical accuracy certification"
  )
  for (bad_matrix_rhs in list(
    cbind(c(0, 0), no_scale_bad_rhs),
    cbind(no_scale_bad_rhs, c(0, 0))
  )) {
    expect_error(
      saddlepoint:::.K2_factor_solve(
        no_scale_factor,
        numeric(2),
        numeric(),
        bad_matrix_rhs
      ),
      "numerical accuracy certification"
    )
  }
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      no_scale_factor,
      numeric(2),
      numeric(),
      matrix(0, 2, 2)
    ),
    matrix(0, 2, 2),
    tolerance = 0
  )

  masked_root <- matrix(0, 3, 3)
  masked_root[1:2, 1:2] <- unstable_root
  masked_root[3, 3] <- 1
  masked_K2 <- tcrossprod(masked_root)
  masked_rhs <- as.vector(masked_K2 %*% c(1, -0.7, 1e30))
  masked_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% masked_root, rep(1, 3))
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      masked_factor,
      numeric(3),
      numeric(),
      masked_rhs
    ),
    "numerical accuracy certification"
  )

  # The terminal decision and solve certificate must not depend on output
  # coordinate order when a tiny covariance coupling is multiplied by a very
  # large right-hand-side component.
  tiny_coupling <- 1e-20
  coupled_correlation <- matrix(c(
    1, 0.3, tiny_coupling,
    0.3, 1, tiny_coupling,
    tiny_coupling, tiny_coupling, 1
  ), 3, byrow = TRUE)
  coupled_root <- diag(c(1e-11, 1e11, 1)) %*%
    t(chol(coupled_correlation))
  coupled_target <- c(1, -0.7, 1e30)
  coordinate_orders <- list(
    c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)
  )
  for (coordinate_order in coordinate_orders) {
    ordered_root <- coupled_root[coordinate_order, , drop = FALSE]
    ordered_K2 <- tcrossprod(ordered_root)
    ordered_rhs <- as.vector(
      ordered_K2 %*% coupled_target[coordinate_order]
    )
    ordered_factor <- local({
      root <- ordered_root
      function(tvec, parameter_vector, A) {
        saddlepoint:::.K2_factor_term(A %*% root, rep(1, 3))
      }
    })
    expect_equal(
      saddlepoint:::.K2_factor_logdet(
        ordered_factor, numeric(3), numeric()
      ),
      as.numeric(determinant(ordered_K2, logarithm = TRUE)$modulus),
      tolerance = 1e-12
    )
    ordinary_target <- (c(1, -0.7, 2) / c(1e-11, 1e11, 1))[
      coordinate_order
    ]
    expect_equal(
      saddlepoint:::.K2_factor_solve(
        ordered_factor, numeric(3), numeric(),
        as.vector(ordered_K2 %*% ordinary_target)
      ),
      ordinary_target,
      tolerance = 1e-12
    )
    expect_error(
      saddlepoint:::.K2_factor_solve(
        ordered_factor, numeric(3), numeric(), ordered_rhs
      ),
      "numerical accuracy certification"
    )
  }

  # Coupled cancellation to an ordinary-scale zero is not itself numerical
  # loss; the mixed certificate is applied after covariance equilibration.
  benign_root <- t(chol(correlation))
  benign_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% benign_root, c(1, 1))
  }
  benign_rhs <- as.vector(correlation %*% c(0, 1))
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      benign_factor, numeric(2), numeric(), benign_rhs
    ),
    c(0, 1),
    tolerance = 1e-14
  )

  ordinary_dimension <- 70L
  ordinary_correlation <- 0.3^abs(outer(
    seq_len(ordinary_dimension), seq_len(ordinary_dimension), `-`
  ))
  ordinary_root <- t(chol(ordinary_correlation))
  ordinary_target <- 1e4 * sin(seq_len(ordinary_dimension))
  ordinary_target[1L] <- 0
  ordinary_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% ordinary_root, rep(1, ordinary_dimension)
    )
  }
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      ordinary_factor, numeric(ordinary_dimension), numeric(),
      as.vector(ordinary_correlation %*% ordinary_target)
    ),
    ordinary_target,
    tolerance = 1e-10
  )

  # Row equilibration must not hide loss introduced when the standardized
  # solution is converted back to the public coordinate scale.
  scaled_dimension <- 5L
  scaled_correlation <- 0.5^abs(outer(
    seq_len(scaled_dimension), seq_len(scaled_dimension), `-`
  ))
  scaled_rows <- 10^seq(-6, 6, length.out = scaled_dimension)
  scaled_root <- scaled_rows * t(chol(scaled_correlation))
  scaled_K2 <- tcrossprod(scaled_root)
  scaled_target <- c(
    0.387119318307123, -1.28520154753104, 0,
    0.188225815574542, 1.26945074078246
  )
  exposed_loss_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% scaled_root, rep(1, scaled_dimension)
    )
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      exposed_loss_factor, numeric(scaled_dimension), numeric(),
      as.vector(scaled_K2 %*% scaled_target)
    ),
    "numerical accuracy certification"
  )

  # Weighted factors can expose a solve error just below the nominal
  # half-precision boundary.  Keep a safety margin in the certificate rather
  # than returning a plausible but inaccurate finite result.
  weighted_root <- matrix(c(
    0x1.ffa27f4c42e1cp-2, -0x1.ecef68e21f7bbp-2,
    0x1.77ad2d90295f1p-1, -0x1.266a450e7ac7ap-1,
    -0x1.65479046f907dp-4, 0x1.f395c25d12d57p-2,
    0x1.bfe68881750cdp+0, -0x1.80a0f56cfee96p+0,
    0x1.989e30b095c24p-2
  ), 3, 3)
  weighted_d <- c(
    0x1.5825fc1eab9a5p+1,
    0x1.c3123764888fcp+17,
    0x1.5b90e5ef373e2p+28
  )
  weighted_rhs <- c(
    0x1.1d57406c5b4fdp+0,
    -0x1.a5dd33fbc225dp-3,
    0x1.ac260169f52a6p+0
  )
  weighted_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% weighted_root, weighted_d)
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      weighted_factor, numeric(3), numeric(), weighted_rhs
    ),
    "numerical accuracy"
  )

  # Exact finite counterexamples must fail closed instead of returning a
  # plausible but badly damaged solve or log-determinant.
  damaged_B <- matrix(c(
    -0x1.7b7ae9113e637p-6, -0x1.c61ea704d14cap+51,
    -0x1.391d3b9ba97b2p-36, 0x1.5254886725324p+20,
    0x1.f365f36d6ce8bp-2, -0x1.9f9a74c87a6b1p+53,
    -0x1.0a412544e31edp-38, -0x1.3e132801de333p+14
  ), 2, 4)
  damaged_d <- c(
    0x1.58b86ef53a2fp-58, 0x1.1721fa532392p+3,
    0x1.488eae360c6c9p-63, 0x1.3a92e82b7a477p+13
  )
  damaged_rhs <- c(
    -0x1.93b54de1c50cbp-47, -0x1.aeaf43e916423p+23
  )
  damaged_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% damaged_B, damaged_d)
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      damaged_factor, numeric(2), numeric(), damaged_rhs
    ),
    "numerical accuracy"
  )

  exact_correlation <- matrix(c(1, 0.5, 0.5, 1), 2)
  exact_root <- t(chol(exact_correlation))
  exact_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% exact_root, c(1, 1))
  }
  expect_error(
    saddlepoint:::.K2_factor_solve(
      exact_factor, numeric(2), numeric(), c(2^53, 2^54)
    ),
    "numerical accuracy"
  )

  damaged_logdet_B <- matrix(c(
    -0x1.ef6921498dc66p-228, 0x1.2d49e2436f5e1p-246,
    0x1.64aadce66705fp-253, -0x1.b1d1d196185f3p-272
  ), 2, 2)
  damaged_logdet_d <- c(
    0x1.72b593fc0b821p+299, 0x1.4391ae14a9428p+344
  )
  damaged_logdet_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% damaged_logdet_B, damaged_logdet_d
    )
  }
  expect_error(
    saddlepoint:::.K2_factor_logdet(
      damaged_logdet_factor, numeric(2), numeric()
    ),
    "log-determinant.*numerical accuracy"
  )

  # Acceptance is a property of the represented covariance, not the raw
  # number of factor columns.  Splitting columns or appending zero columns
  # therefore leaves both terminal results unchanged.
  split_count <- 10000L
  split_B <- exact_root[, rep(1:2, each = split_count), drop = FALSE]
  split_d <- rep(1 / split_count, 2L * split_count)
  split_B <- cbind(split_B, matrix(0, 2, 10))
  split_d <- c(split_d, rep(1, 10))
  split_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(A %*% split_B, split_d)
  }
  expect_equal(
    saddlepoint:::.K2_factor_logdet(
      split_factor, numeric(2), numeric()
    ),
    log(0.75),
    tolerance = 1e-11
  )
  expect_equal(
    saddlepoint:::.K2_factor_solve(
      split_factor, numeric(2), numeric(), c(0.2, -0.3)
    ),
    solve(exact_correlation, c(0.2, -0.3)),
    tolerance = 1e-11
  )

  # Two Gram products can share the same dot-product rounding loss.  Certify
  # factor formation itself so a split, ill-conditioned factor cannot pass
  # merely because its rounded covariance is self-consistent.
  split_condition <- 5100000
  split_correlation <-
    (split_condition - 1) / (split_condition + 1)
  split_K2 <- matrix(c(
    1, split_correlation,
    split_correlation, 1
  ), 2)
  split_root <- t(chol(split_K2))
  replicated_columns <- 500L
  rounded_B <- split_root[
    , rep(1:2, each = replicated_columns), drop = FALSE
  ] / sqrt(replicated_columns)
  rounded_factor <- function(tvec, parameter_vector, A) {
    saddlepoint:::.K2_factor_term(
      A %*% rounded_B, rep(1, 2L * replicated_columns)
    )
  }
  rounded_rhs <- as.vector(split_K2 %*% c(1, -0.7))
  expect_error(
    saddlepoint:::.K2_factor_solve(
      rounded_factor, numeric(2), numeric(), rounded_rhs
    ),
    "numerical accuracy"
  )
  expect_error(
    saddlepoint:::.K2_factor_logdet(
      rounded_factor, numeric(2), numeric()
    ),
    "numerical accuracy"
  )

  # A structurally exact zero solution coordinate remains certifiable.
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

  theta_multinomial <- c(1, 0.2, 0.3, 0.5)
  multinomial_K2 <- as.matrix(
    MultinomialCGF$K2(rep(0, 3), theta_multinomial)
  )
  eig <- eigen(multinomial_K2, symmetric = TRUE)
  positive <- eig$values > 100 * .Machine$double.eps
  multinomial_root <- eig$vectors[, positive, drop = FALSE] %*%
    diag(sqrt(eig$values[positive]), sum(positive))
  map_matrix <- unstable_root %*%
    solve(crossprod(multinomial_root), t(multinomial_root))
  unstable_mapped_multinomial <- linearlyMappedCGF(
    MultinomialCGF,
    map_matrix,
    iidReps = 1L
  )
  mapped_K2 <- as.matrix(
    unstable_mapped_multinomial$K2(c(0, 0), theta_multinomial)
  )
  mapped_rhs <- as.vector(mapped_K2 %*% c(1, -0.7))
  expect_error(
    unstable_mapped_multinomial$K2_solve(
      c(0, 0),
      theta_multinomial,
      mapped_rhs
    ),
    "numerical accuracy certification"
  )

  # Completing the model with an independent stable coordinate and then
  # applying an identity map must not bypass the same certification failure.
  fixed_poisson <- PoissonModelCGF(
    lambda = function(p) 1 + 0 * sum(p),
    iidReps = 1L
  )
  completed <- concatenationCGF(
    list(unstable_mapped_multinomial, fixed_poisson),
    component_dims = c(2L, 1L),
    iidReps = 1L
  )
  identity_mapped <- linearlyMappedCGF(
    completed,
    diag(3),
    iidReps = 1L
  )
  completed_tvec <- numeric(3)
  completed_rhs <- as.vector(
    identity_mapped$K2(completed_tvec, theta_multinomial) %*%
      c(1, -0.7, 1e30)
  )
  for (model in list(completed, identity_mapped)) {
    expect_error(
      model$K2_solve(
        completed_tvec,
        theta_multinomial,
        completed_rhs
      ),
      "numerical accuracy certification"
    )
  }

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

test_that("solve certification refines conservative condition metadata", {
  # The determinant/trace proof is intentionally conservative.  It is enough
  # to certify a log-determinant, but an accurate row-scaled solve can require
  # the exact SVD metadata in the componentwise solution certificate.
  m <- 7L
  correlation <- 0x1.7a917fb1e6666p-1
  R <- chol(correlation^abs(outer(seq_len(m), seq_len(m), `-`)))
  scale <- c(
    0x1.6e57fae5a3805p-10, 0x1.0064108ec234p+14,
    0x1.d73f5e463ded1p-19, 0x1.5f827a39ab9fep-4,
    0x1.087b56d22de4dp+39, 0x1.b863eb1a57085p-19,
    0x1.3fd3da381a1e7p+1
  )
  z <- c(0, 0, 0, 0x1.2d01c7de1d8d2p-2, 0, 0, 0)
  discrepancy <- rep(0x1.85f0318b81401p-57, m)
  w <- solve(R, solve(t(R), z))
  solution <- w / scale

  terminal <- saddlepoint:::.K2_terminal_root_guard(
    R, discrepancy, covariance_is_ad = FALSE
  )
  expect_equal(as.numeric(terminal$condition_is_bound), 1)
  expect_equal(saddlepoint:::.K2_solution_guard(
    solution, R, w, z, scale, discrepancy,
    terminal$rho, terminal$lambda_min, terminal$condition_is_bound
  ), solution, tolerance = 0)
})

test_that("scale-aware solve rejection is reentrant on AD tapes", {
  correlation_root <- t(chol(matrix(c(1, 0.3, 0.3, 1), 2L)))
  factor_root <- function(p) {
    diag(c(10^-p[1], 10^p[1])) %*% correlation_root
  }
  cgf <- createCGF(
    K = function(tvec, p) 0 * p[1],
    K1 = function(tvec, p) 0 * tvec,
    K2 = function(tvec, p) tcrossprod(factor_root(p)),
    K3operator = function(tvec, p, a, b, c) 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) 0 * p[1],
    K2_factor = function(tvec, p, A) {
      saddlepoint:::.K2_factor_term(
        A %*% factor_root(p),
        c(1, 1)
      )
    },
    K2_factor_terminal = function() TRUE
  )
  target <- cbind(c(1, -0.7), c(-0.2, 0.4))
  objective <- function(p) {
    K2 <- cgf$K2(numeric(2), p)
    rhs <- K2 %*% target
    sum(cgf$K2_solve(numeric(2), p, rhs))
  }
  tape <- RTMB::MakeTape(objective, 0)
  derivative_tape <- tape$jacfun()
  evaluate <- function(p) {
    c(
      value = tape(p),
      gradient = tape$jacobian(p),
      hessian = derivative_tape$jacobian(p)
    )
  }

  expected <- c(value = 0.5, gradient = 0, hessian = 0)
  expect_equal(evaluate(0), expected, tolerance = 1e-12)
  expect_true(all(is.nan(evaluate(11))))
  expect_equal(evaluate(0), expected, tolerance = 1e-12)
})

test_that("factor solve rejects unreliable Hessians and recovers", {
  covariance_root <- function(p) {
    zero <- 0 * p[1]
    root <- diag(2) + zero
    root[1, 2] <- zero
    root[2, 1] <- 1 - p[1]
    root[2, 2] <- sqrt(p[1] * (2 - p[1]))
    root
  }
  factor <- function(tvec, p, A) {
    saddlepoint:::.K2_factor_term(
      A %*% covariance_root(p), c(1, 1)
    )
  }
  rhs <- c(0.3, -0.2)
  objective <- function(p) {
    sum(saddlepoint:::.K2_factor_solve(
      factor, numeric(2), p, rhs
    ))
  }

  # The value-only calculation remains accurate near the boundary.  The same
  # covariance cannot safely support a generic second derivative in double
  # precision and must invalidate the whole taped evaluation.
  expect_equal(objective(1e-5), 0.1 / (2 - 1e-5), tolerance = 1e-10)
  tape <- RTMB::MakeTape(objective, 0.1)
  derivative_tape <- tape$jacfun()
  evaluate <- function(p) {
    c(
      value = tape(p),
      gradient = tape$jacobian(p),
      hessian = derivative_tape$jacobian(p)
    )
  }
  expected <- c(
    value = 0.1 / 1.9,
    gradient = 0.1 / 1.9^2,
    hessian = 0.2 / 1.9^3
  )
  expect_equal(evaluate(0.1), expected, tolerance = 1e-10)
  moderately_conditioned_expected <- c(
    value = 0.1 / (2 - 0.003),
    gradient = 0.1 / (2 - 0.003)^2,
    hessian = 0.2 / (2 - 0.003)^3
  )
  expect_equal(
    evaluate(0.003), moderately_conditioned_expected,
    tolerance = 1e-7
  )

  fixed_condition <- 9000
  fixed_correlation <-
    (fixed_condition - 1) / (fixed_condition + 1)
  fixed_K2 <- matrix(c(
    1, fixed_correlation,
    fixed_correlation, 1
  ), 2)
  fixed_root <- t(chol(fixed_K2))
  fixed_factor <- function(tvec, p, A) {
    saddlepoint:::.K2_factor_term(A %*% fixed_root, c(1, 1))
  }
  fixed_objective <- function(p) {
    target <- c(p[1], p[1]^2)
    sum(saddlepoint:::.K2_factor_solve(
      fixed_factor, numeric(2), p,
      as.vector(fixed_K2 %*% target)
    ))
  }
  fixed_tape <- RTMB::MakeTape(fixed_objective, 0.2)
  expect_equal(
    c(
      fixed_tape(0.2),
      fixed_tape$jacobian(0.2),
      fixed_tape$jacfun()$jacobian(0.2)
    ),
    c(0.24, 1.4, 2),
    tolerance = 1e-8
  )
  expect_true(all(is.nan(evaluate(1e-5))))
  expect_equal(evaluate(0.1), expected, tolerance = 1e-10)

  flat_factor <- function(tvec, p, A) {
    saddlepoint:::.K2_factor_term(
      A %*% covariance_root(1e-5 + p[1]^2), c(1, 1)
    )
  }
  flat_objective <- function(p) {
    sum(saddlepoint:::.K2_factor_solve(
      flat_factor, numeric(2), p, rhs
    ))
  }
  flat_tape <- RTMB::MakeTape(flat_objective, 0.1)
  flat_derivative_tape <- flat_tape$jacfun()
  flat_evaluate <- function(p) {
    c(
      flat_tape(p),
      flat_tape$jacobian(p),
      flat_derivative_tape$jacobian(p)
    )
  }
  expect_true(all(is.nan(flat_evaluate(0))))
  expect_true(all(is.finite(flat_evaluate(0.1))))
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
  # but does not certify its missing counterpart for an enclosing RSS update.
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
