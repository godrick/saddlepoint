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
