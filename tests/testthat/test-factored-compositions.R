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
    func_T = function(tvec, p) 5 + p[1]^3
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
  expect_equal(combined$.private_api$func_T(tvec, theta), 0)
  overridden <- sumOfIndependentCGF(
    list(child),
    iidReps = 1L,
    func_T = function(tvec, p) 77 + 0 * p[1]
  )
  expect_equal(overridden$.private_api$func_T(tvec, theta), 77)
})

test_that("independent sums reject cancellation and recover on one tape", {
  make_mapped_poisson <- function(coefficient) {
    linearlyMappedCGF(
      PoissonModelCGF(lambda = function(p) 1 + 0 * p[1], iidReps = 1L),
      matrix(coefficient, 1L, 1L),
      iidReps = 1L
    )
  }
  lossy <- sumOfIndependentCGF(
    lapply(c(2^18, 1, -2^18), make_mapped_poisson),
    iidReps = 1L
  )
  tape <- RTMB::MakeTape(function(p) {
    lossy$K3operator(0, p, p[1], p[1], p[1])
  }, 0)

  expect_true(all(is.nan(c(
    tape(1),
    tape$jacobian(1),
    tape$jacfun()$jacobian(1)
  ))))
  expect_equal(c(
    tape(0),
    tape$jacobian(0),
    tape$jacfun()$jacobian(0)
  ), c(0, 0, 0))
})
