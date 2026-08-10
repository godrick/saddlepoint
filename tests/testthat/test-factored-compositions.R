make_factored_composition_child <- function() {
  cumulant <- function(tvec, parameter_vector) {
    exp(parameter_vector[1] + tvec)
  }

  createCGF(
    K = function(tvec, parameter_vector) {
      sum(exp(parameter_vector[1]) * (exp(tvec) - 1))
    },
    K1 = function(tvec, parameter_vector) cumulant(tvec, parameter_vector),
    K2 = function(tvec, parameter_vector) {
      diag(cumulant(tvec, parameter_vector), length(tvec))
    },
    K3operator = function(tvec, parameter_vector, a, b, c) {
      sum(cumulant(tvec, parameter_vector) * a * b * c)
    },
    K4operator = function(tvec, parameter_vector, a, b, c, d) {
      sum(cumulant(tvec, parameter_vector) * a * b * c * d)
    },
    op_name = "factoredCompositionChild"
  )
}

factored_tape_vgh <- function(objective, theta) {
  tape <- RTMB::MakeTape(objective, theta)
  c(
    value = tape(theta),
    gradient = tape$jacobian(theta),
    hessian = as.vector(tape$jacfun()$jacobian(theta))
  )
}

expect_factored_composition_vgh <- function(cgf, tvec, theta, factor_A,
                                             factor_d) {
  dense_Q <- function(p) {
    A <- factor_A(p)
    A %*% (factor_d(p) * t(A))
  }

  for (methods in list(
    c("K3K3operatorAABBCC_factored", "K3K3operatorAABBCC"),
    c("K3K3operatorABCABC_factored", "K3K3operatorABCABC")
  )) {
    candidate <- RTMB::MakeTape(function(p) {
      cgf$.private_api[[methods[[1L]]]](
        tvec, p, factor_A(p), factor_d(p)
      )
    }, theta)
    reference <- RTMB::MakeTape(function(p) {
      cgf[[methods[[2L]]]](tvec, p, dense_Q(p))
    }, theta)

    expect_equal(candidate(theta), reference(theta), tolerance = 1e-10)
    expect_equal(
      candidate$jacobian(theta),
      reference$jacobian(theta),
      tolerance = 1e-9
    )
    expect_equal(
      candidate$jacfun()$jacobian(theta),
      reference$jacfun()$jacobian(theta),
      tolerance = 1e-8
    )
  }
}

test_that("multi-block IID preserves thin contractions and zero weights", {
  block_dim <- 4L
  n_blocks <- 3L
  total_dim <- block_dim * n_blocks
  theta <- c(0.2, 0)
  cgf <- iidReplicatesCGF(
    make_factored_composition_child(),
    iidReps = n_blocks,
    block_size = block_dim
  )
  tvec <- seq(-0.08, 0.09, length.out = total_dim)

  set.seed(941)
  A0 <- matrix(rnorm(2L * total_dim), total_dim, 2L)
  A_delta <- matrix(
    seq(-0.03, 0.03, length.out = 2L * total_dim),
    total_dim,
    2L
  )
  factor_A <- function(p) A0 + p[2] * A_delta
  factor_d <- function(p) c(exp(0.1 * p[1]), p[2]^2)

  expect_factored_composition_vgh(cgf, tvec, theta, factor_A, factor_d)

  set.seed(951)
  map <- matrix(rnorm(2L * total_dim), nrow = 2L) / sqrt(total_dim)
  mapped_iid <- linearlyMappedCGF(cgf, map, iidReps = 1L)
  mapped_flat <- linearlyMappedCGF(
    make_factored_composition_child(), map, iidReps = 1L
  )
  mapped_tvec <- c(0.02, -0.03)
  iid_tape <- RTMB::MakeTape(
    function(p) mapped_iid$.private_api$func_T(mapped_tvec, p), theta
  )
  flat_tape <- RTMB::MakeTape(
    function(p) mapped_flat$.private_api$func_T(mapped_tvec, p), theta
  )
  expect_equal(
    c(
      iid_tape(theta),
      iid_tape$jacobian(theta),
      iid_tape$jacfun()$jacobian(theta)
    ),
    c(
      flat_tape(theta),
      flat_tape$jacobian(theta),
      flat_tape$jacfun()$jacobian(theta)
    ),
    tolerance = 1e-9
  )
})

test_that("concatenation preserves thin contractions and zero weights", {
  n_components <- 3L
  component_dim <- 4L
  total_dim <- n_components * component_dim
  theta <- c(0.13, 0)
  child <- make_factored_composition_child()
  cgf <- concatenationCGF(
    rep(list(child), n_components),
    component_dims = component_dim,
    iidReps = 1L
  )
  tvec <- seq(-0.04, 0.05, length.out = total_dim)

  set.seed(7301)
  A0 <- matrix(rnorm(2L * total_dim), total_dim, 2L) / sqrt(total_dim)
  A_delta <- matrix(
    seq(-0.03, 0.03, length.out = 2L * total_dim),
    total_dim,
    2L
  )
  factor_A <- function(p) A0 + p[2] * A_delta
  factor_d <- function(p) c(exp(0.1 * p[1]), p[2]^2)

  expect_factored_composition_vgh(cgf, tvec, theta, factor_A, factor_d)
})

test_that("singleton wrappers preserve authoritative correction methods", {
  child <- createCGF(
    K = function(tvec, p) 0.5 * exp(p[1]) * sum(tvec * tvec),
    K1 = function(tvec, p) exp(p[1]) * tvec,
    K2 = function(tvec, p) diag(exp(p[1]), length(tvec)),
    K3operator = function(tvec, p, a, b, c) 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) 0 * p[1],
    K4operatorAABB_factored = function(tvec, p, A, d) 8 + p[1]^2,
    K3K3operatorAABBCC_factored = function(tvec, p, A, d) 8 - p[1],
    K3K3operatorABCABC_factored = function(tvec, p, A, d) 12 + p[1],
    func_T = function(tvec, p) 5 + p[1]^3,
    neg_ll = function(tvec, p) 41 + p[1]^2,
    K2_solve = function(tvec, p, rhs) rhs / (3 + p[1]),
    logdetK2 = function(tvec, p) length(tvec) * log(3 + p[1]),
    analytic_tvec_hat = function(x, p) x / (3 + p[1])
  )
  singleton <- sumOfIndependentCGF(list(child), iidReps = 1L)
  tvec <- 0.2
  theta <- 0.4
  A <- matrix(1, 1L, 1L)

  for (method_name in c(
    "K4operatorAABB_factored",
    "K3K3operatorAABBCC_factored",
    "K3K3operatorABCABC_factored"
  )) {
    expect_equal(
      singleton$.private_api[[method_name]](tvec, theta, A, 1),
      child$.private_api[[method_name]](tvec, theta, A, 1),
      tolerance = 0
    )
  }

  child_tape <- RTMB::MakeTape(
    function(p) child$.private_api$func_T(tvec, p), theta
  )
  singleton_tape <- RTMB::MakeTape(
    function(p) singleton$.private_api$func_T(tvec, p), theta
  )
  expect_equal(
    c(
      singleton_tape(theta),
      singleton_tape$jacobian(theta),
      singleton_tape$jacfun()$jacobian(theta)
    ),
    c(
      child_tape(theta),
      child_tape$jacobian(theta),
      child_tape$jacfun()$jacobian(theta)
    ),
    tolerance = 1e-12
  )
  expect_equal(
    singleton$K2_solve(tvec, theta, 2),
    child$K2_solve(tvec, theta, 2),
    tolerance = 0
  )
  expect_equal(singleton$logdetK2(tvec, theta), child$logdetK2(tvec, theta))
  expect_equal(
    singleton$.private_api$neg_ll(tvec, theta),
    child$.private_api$neg_ll(tvec, theta)
  )
  expect_true(singleton$has_analytic_tvec_hat)
  expect_equal(
    singleton$analytic_tvec_hat(2, theta),
    child$analytic_tvec_hat(2, theta)
  )

  stages <- list(
    adaptCGF(child, function(p) p),
    .exponentialTiltCGF_internal(child, function(p) 0 * p[1]),
    shiftedCGF(child, 0),
    sumOfiidCGF(child, n = 1),
    linearlyMappedCGF(child, matrix(1, 1L, 1L), iidReps = 1L)
  )
  for (stage in stages) {
    expect_true(all(vapply(
      c(
        "K4operatorAABB_factored",
        "K3K3operatorAABBCC_factored",
        "K3K3operatorABCABC_factored"
      ),
      function(name) saddlepoint:::.factored_delegate_is_safe(
        stage$.private_api[[name]]
      ),
      logical(1)
    )))
  }

  combined <- sumOfIndependentCGF(list(child, child), iidReps = 1L)
  expect_equal(
    combined$.private_api$func_T(tvec, theta),
    2 * (8 + theta^2) / 8,
    tolerance = 1e-12
  )
  overridden <- sumOfIndependentCGF(
    list(child),
    iidReps = 1L,
    func_T = function(tvec, p) 77 + 0 * p[1]
  )
  expect_equal(overridden$.private_api$func_T(tvec, theta), 77)

  numerical_override <- sumOfIndependentCGF(
    list(child),
    iidReps = 1L,
    K2_solve = function(tvec, p, rhs) 7 * rhs + 0 * p[1],
    logdetK2 = function(tvec, p) 13 + 0 * p[1],
    neg_ll = function(tvec, p) 17 + 0 * p[1],
    analytic_tvec_hat = function(x, p) x + 19 + 0 * p[1]
  )
  expect_equal(numerical_override$K2_solve(tvec, theta, 2), 14)
  expect_equal(numerical_override$logdetK2(tvec, theta), 13)
  expect_equal(numerical_override$.private_api$neg_ll(tvec, theta), 17)
  expect_equal(numerical_override$analytic_tvec_hat(2, theta), 21)
})

test_that("independent sums add child factored K4 methods directly", {
  calls <- new.env(parent = emptyenv())
  calls$primitive <- 0L
  calls$factored <- 0L
  make_child <- function() createCGF(
    K = function(tvec, p) 0 * p[1],
    K1 = function(tvec, p) numeric(length(tvec)) * p[1],
    K2 = function(tvec, p) diag(length(tvec)) * (1 + 0 * p[1]),
    K3operator = function(tvec, p, a, b, c) 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) {
      calls$primitive <- calls$primitive + 1L
      0 * p[1]
    },
    K4operatorAABB_factored = function(tvec, p, A, d) {
      calls$factored <- calls$factored + 1L
      8 + p[1]^2
    }
  )
  cgf <- sumOfIndependentCGF(
    list(make_child(), make_child()), iidReps = 1L
  )
  tvec <- c(0.1, -0.2, 0.3)
  A <- matrix(c(1, 0.2, -0.1, 0.5, 0.8, 0.3), nrow = 3L)
  d <- c(1, 0.7)

  tape <- RTMB::MakeTape(
    function(p) cgf$.private_api$K4operatorAABB_factored(tvec, p, A, d),
    0.4
  )
  expect_equal(
    c(tape(0.4), tape$jacobian(0.4), tape$jacfun()$jacobian(0.4)),
    c(16.32, 1.6, 4),
    tolerance = 1e-12
  )
  expect_identical(calls$primitive, 0L)
  expect_identical(calls$factored, 2L)
})

test_that("independent sums keep ordinary near-cancellation differentiable", {
  make_mapped_poisson <- function(coefficient) {
    linearlyMappedCGF(
      PoissonModelCGF(lambda = function(p) exp(p[1]), iidReps = 1L),
      matrix(coefficient, 1L, 1L),
      iidReps = 1L
    )
  }
  coefficients <- c(1, -(1 - 2^-52))
  nearly_symmetric <- sumOfIndependentCGF(
    lapply(coefficients, make_mapped_poisson),
    iidReps = 1L
  )
  tape <- RTMB::MakeTape(function(p) {
    nearly_symmetric$K3operator(0, p, 1, 1, 1)
  }, 0.2)

  expected <- sum(coefficients^3) * exp(0.2)
  result <- c(tape(0.2), tape$jacobian(0.2), tape$jacfun()$jacobian(0.2))
  expect_true(all(is.finite(result)))
  expect_equal(result, rep(expected, 3), tolerance = 1e-14)
})

test_that("K4 delegation preserves equivalent zero-column factors", {
  theta <- c(0.1, 0)
  tvec <- c(-0.02, 0.005, 0.03)
  v <- c(
    0.213962502184879833,
    0.47965813457087475,
    0.087828704973743787
  )
  factor_zero <- cbind(diag(3), numeric(3))
  factor_delta <- cbind(matrix(0, 3, 3), v)
  column_factor <- function(p) factor_zero + p[2] * factor_delta
  column_weight <- function(p) rep(1, 4)
  weight_factor <- function(p) cbind(diag(3), v)
  weight_weight <- function(p) {
    c(1, 1, 1, 0) + c(0, 0, 0, 1) * p[2]^2
  }
  dense_Q <- function(p) diag(3) + p[2]^2 * tcrossprod(v)

  independent_sum <- sumOfIndependentCGF(
    list(
      make_factored_composition_child(),
      make_factored_composition_child()
    ),
    iidReps = 1L
  )
  sum_column <- factored_tape_vgh(function(p) {
    independent_sum$.private_api$K4operatorAABB_factored(
      tvec, p, column_factor(p), column_weight(p)
    )
  }, theta)
  sum_weight <- factored_tape_vgh(function(p) {
    independent_sum$.private_api$K4operatorAABB_factored(
      tvec, p, weight_factor(p), weight_weight(p)
    )
  }, theta)
  expect_equal(sum_column, sum_weight, tolerance = 1e-12)
  expect_equal(
    sum_column[["hessian4"]],
    8 * sum(exp(theta[1] + tvec) * v^2),
    tolerance = 1e-12
  )

  count <- PoissonModelCGF(
    lambda = function(p) exp(p[1]),
    iidReps = 1L
  )
  rss <- randomlyStoppedSumCGF(
    count,
    make_factored_composition_child(),
    block_size = 3L,
    iidReps = 1L
  )
  rss_column <- factored_tape_vgh(function(p) {
    rss$.private_api$K4operatorAABB_factored(
      tvec, p, column_factor(p), column_weight(p)
    )
  }, theta)
  rss_weight <- factored_tape_vgh(function(p) {
    rss$.private_api$K4operatorAABB_factored(
      tvec, p, weight_factor(p), weight_weight(p)
    )
  }, theta)
  rss_dense <- factored_tape_vgh(function(p) {
    rss$K4operatorAABB(tvec, p, dense_Q(p))
  }, theta)
  expect_equal(rss_column, rss_weight, tolerance = 1e-11)
  expect_equal(rss_column, rss_dense, tolerance = 1e-11)
})

test_that("mapped factored bridges accept fixed Matrix factors with AD weights", {
  dense_child <- createCGF(
    K = function(tvec, p) 0 * p[1],
    K1 = function(tvec, p) 0 * tvec + 0 * p[1],
    K2 = function(tvec, p) diag(length(tvec)) + 0 * p[1],
    K3operator = function(tvec, p, a, b, c) sum(a * b * c) + 0 * p[1],
    K4operator = function(tvec, p, a, b, c, d) {
      sum(a * b * c * d) + 0 * p[1]
    },
    K4operatorAABB = function(tvec, p, Q) sum(Q^2) + p[1],
    K3K3operatorAABBCC = function(tvec, p, Q) sum(Q)^2 + p[1],
    K3K3operatorABCABC = function(tvec, p, Q) sum(Q^3) + p[1]
  )
  mapped <- linearlyMappedCGF(dense_child, diag(3), iidReps = 1L)
  factor_A <- Matrix::Matrix(
    cbind(diag(3), c(0.2, 0.4, 0.1)),
    sparse = FALSE
  )
  theta <- c(0.1, 0)
  factor_d <- function(p) {
    c(1, 1, 1, 0) + c(0, 0, 0, 1) * p[2]^2
  }

  for (method_name in c(
    "K4operatorAABB_factored",
    "K3K3operatorAABBCC_factored",
    "K3K3operatorABCABC_factored"
  )) {
    mapped_vgh <- factored_tape_vgh(function(p) {
      mapped$.private_api[[method_name]](
        numeric(3), p, factor_A, factor_d(p)
      )
    }, theta)
    child_vgh <- factored_tape_vgh(function(p) {
      dense_child$.private_api[[method_name]](
        numeric(3), p, base::as.matrix(factor_A), factor_d(p)
      )
    }, theta)
    expect_equal(mapped_vgh, child_vgh, tolerance = 1e-12)
  }
})

test_that("mapped dense contractions preserve AD with numeric Matrix inputs", {
  child <- PoissonModelCGF(
    lambda = adaptor(indices = 1L),
    iidReps = 2L
  )
  A <- matrix(c(
    1, 0.2,
    0.3, -0.1,
    0.8, 0.4
  ), nrow = 3L, byrow = TRUE)
  mapped <- linearlyMappedCGF(child, A)
  tvec <- c(0.01, -0.02, 0.03)
  Q <- matrix(c(
    2, 0.2, 0.1,
    0.2, 1.5, 0.3,
    0.1, 0.3, 1.2
  ), 3L, 3L)
  theta <- 1.1
  inner_tvec <- as.vector(t(A) %*% tvec)
  inner_Q <- t(A) %*% Q %*% A

  for (Q_input in list(Q, Matrix::Matrix(Q, sparse = FALSE))) {
    for (method in c("K3K3operatorAABBCC", "K3K3operatorABCABC")) {
      candidate <- factored_tape_vgh(
        function(p) mapped[[method]](tvec, p, Q_input), theta
      )
      reference <- factored_tape_vgh(
        function(p) child[[method]](inner_tvec, p, inner_Q), theta
      )
      expect_equal(candidate, reference, tolerance = 1e-12)
    }
  }
})

test_that("parameter-dependent maps accept numeric Matrix contractions", {
  child <- PoissonModelCGF(
    lambda = adaptor(indices = 1L),
    iidReps = 2L
  )
  A0 <- matrix(c(1, 0.2, 0.3, -0.1, 0.8, 0.4), 3L, 2L)
  A_delta <- matrix(c(0.02, -0.01, 0.03, 0.01, -0.02, 0.01), 3L, 2L)
  A_fun <- function(p) A0 + p[1] * A_delta
  mapped <- linearlyMappedCGF(child, A_fun)
  tvec <- c(0.01, -0.02, 0.03)
  Q <- matrix(c(2, 0.2, 0.1, 0.2, 1.5, 0.3, 0.1, 0.3, 1.2), 3L)
  Q_matrix <- Matrix::Matrix(Q, sparse = FALSE)
  theta <- 1.1

  for (method in c("K3K3operatorAABBCC", "K3K3operatorABCABC")) {
    candidate <- factored_tape_vgh(
      function(p) mapped[[method]](tvec, p, Q_matrix), theta
    )
    reference <- factored_tape_vgh(function(p) {
      A <- A_fun(p)
      child[[method]](
        as.vector(t(A) %*% tvec), p,
        t(A) %*% Q %*% A
      )
    }, theta)
    expect_equal(candidate, reference, tolerance = 1e-11)
  }

  operator_B <- matrix(c(1, 0.2, 0.3, 0.4, 0.5, 0.6), 2L, 3L)
  expect_equal(
    factored_tape_vgh(
      function(p) sum(mapped$K2operatorAK2AT(
        tvec, p, Matrix::Matrix(operator_B, sparse = FALSE)
      )), theta
    ),
    factored_tape_vgh(
      function(p) sum(mapped$K2operatorAK2AT(tvec, p, operator_B)), theta
    ),
    tolerance = 1e-12
  )

  factor_B <- matrix(c(1, 0.2, 0.3, 0.4, 0.5, 0.6), 3L, 2L)
  factor_d <- function(p) {
    d <- p[rep(1L, 2L)] * 0
    d[1L] <- 1
    d[2L] <- p[1]^2
    d
  }
  for (method in c(
    "K4operatorAABB_factored",
    "K3K3operatorAABBCC_factored",
    "K3K3operatorABCABC_factored"
  )) {
    candidate <- factored_tape_vgh(function(p) {
      mapped$.private_api[[method]](
        tvec, p, Matrix::Matrix(factor_B, sparse = FALSE), factor_d(p)
      )
    }, theta)
    reference <- factored_tape_vgh(function(p) {
      mapped$.private_api[[method]](tvec, p, factor_B, factor_d(p))
    }, theta)
    expect_equal(candidate, reference, tolerance = 1e-11)
  }

  multinomial_map <- linearlyMappedCGF(MultinomialCGF, A_fun)
  mapped_factor <- saddlepoint:::.K2_factor_method(multinomial_map)
  multinomial_parameter <- c(10, 0.4, 0.6)
  factor_value <- function(p, B) {
    terms <- mapped_factor(tvec, p, B)
    sum(terms[[1L]]$B) + sum(terms[[1L]]$d)
  }
  expect_equal(
    factored_tape_vgh(
      function(p) factor_value(
        p, Matrix::Matrix(operator_B, sparse = FALSE)
      ), multinomial_parameter
    ),
    factored_tape_vgh(
      function(p) factor_value(p, operator_B), multinomial_parameter
    ),
    tolerance = 1e-11
  )
})

test_that("mapped factored contractions pass ordinary numeric factors", {
  strict_factored <- function(value) {
    force(value)
    function(tvec, p, B, d) {
      if (!is.matrix(B)) stop("factor B must be an ordinary matrix")
      value + 0 * sum(B) + 0 * sum(d) + 0 * sum(p)
    }
  }
  child <- createCGF(
    K = function(tvec, p) sum(exp(p[1]) * (exp(tvec) - 1)),
    K1 = function(tvec, p) exp(p[1] + tvec),
    K2 = function(tvec, p) diag(exp(p[1] + tvec), length(tvec)),
    K3operator = function(tvec, p, a, b, c) {
      sum(exp(p[1] + tvec) * a * b * c)
    },
    K4operator = function(tvec, p, a, b, c, d) {
      sum(exp(p[1] + tvec) * a * b * c * d)
    },
    K4operatorAABB_factored = strict_factored(8),
    K3K3operatorAABBCC_factored = strict_factored(8),
    K3K3operatorABCABC_factored = strict_factored(12)
  )
  mapped <- linearlyMappedCGF(
    child,
    matrix(c(1, 0.2, 0.3, -0.1, 0.8, 0.4), nrow = 2L, byrow = TRUE)
  )

  expect_equal(mapped$K4operatorAABB(c(0, 0), 1, diag(2)), 8)
  expect_equal(mapped$K3K3operatorAABBCC(c(0, 0), 1, diag(2)), 8)
  expect_equal(mapped$K3K3operatorABCABC(c(0, 0), 1, diag(2)), 12)
})

test_that("shared-Poisson IID correction remains finite", {
  dimension <- 3L
  blocks <- 5L
  shared_scalar <- PoissonModelCGF(
    lambda = adaptor(indices = 2), iidReps = 1L
  )
  shared_vector <- linearlyMappedCGF(
    shared_scalar, matrix(1, nrow = dimension), iidReps = 1L
  )
  independent_vector <- PoissonModelCGF(
    lambda = adaptor(indices = 1), iidReps = dimension
  )
  model <- iidReplicatesCGF(
    sumOfIndependentCGF(
      list(independent_vector, shared_vector), iidReps = 1L
    ),
    iidReps = blocks
  )
  tvec <- seq(-0.04, 0.05, length.out = dimension * blocks)
  theta <- c(14, 7)
  tape <- RTMB::MakeTape(
    function(p) model$.private_api$func_T(tvec, p), theta
  )
  result <- c(
    tape(theta),
    tape$jacobian(theta),
    tape$jacfun()$jacobian(theta)
  )

  expect_true(all(is.finite(result)))
  expect_equal(
    result,
    c(
      -0.049469809893120435,
      0.001190996049301184,
      0.004685123600414837,
      6.723910735767176e-05,
      -4.747628002299675e-04,
      -4.747628002299674e-04,
      -3.890811425157328e-04
    ),
    tolerance = 1e-11
  )
})
