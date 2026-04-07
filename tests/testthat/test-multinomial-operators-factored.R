
make_Q_factored <- function(d, rank = max(1L, d - 1L), seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(d * rank), nrow = d, ncol = rank)
  dd <- runif(rank, min = 0.4, max = 1.2)
  Q <- A %*% diag(dd, nrow = rank) %*% t(A)
  list(Q = Q, A = A, d = dd)
}

expect_multinomial_factored_matches_dense <- function(cgf, param, tvec, seed) {
  fd <- make_Q_factored(length(tvec), rank = 2, seed = seed)

  expect_equal(
    as.numeric(cgf$.private_api$K4operatorAABB_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K4operatorAABB(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorAABBCC_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K3K3operatorAABBCC(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorABCABC_factored(tvec, param, fd$A, fd$d)),
    as.numeric(cgf$K3K3operatorABCABC(tvec, param, fd$Q)),
    tolerance = 1e-10
  )
}

test_that("MultinomialFamilyCGF factored operators match explicit dense-Q operators", {
  cgf <- MultinomialFamilyCGF$new()

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)
  tvec <- c(0.05, -0.08, 0.02)

  expect_multinomial_factored_matches_dense(cgf, param, tvec, seed = 10)
})


test_that("K4operatorAABB_factored matches explicit K4operatorAABB for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)
  tvec <- c(0.05, -0.08, 0.02)

  fd <- make_Q_factored(3, rank = 2, seed = 10)

  K4_fact <- cgf$.private_api$K4operatorAABB_factored

  got_fact <- K4_fact(tvec, param, fd$A, fd$d)
  got_Q    <- cgf$K4operatorAABB(tvec, param, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})


test_that("K3K3operatorAABBCC_factored matches explicit K3K3operatorAABBCC for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 8
  odds <- c(0.2, 0.3, 0.5)
  param <- c(N, odds)
  tvec <- c(0.02, -0.03, 0.01)

  fd <- make_Q_factored(3, rank = 2, seed = 20)

  K3K3_fact <- cgf$.private_api$K3K3operatorAABBCC_factored

  got_fact <- K3K3_fact(tvec, param, fd$A, fd$d)
  got_Q    <- cgf$K3K3operatorAABBCC(tvec, param, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})


test_that("K3K3operatorABCABC_factored matches explicit K3K3operatorABCABC for Q=A diag(d) A^T", {
  cgf <- MultinomialCGF

  N <- 9
  odds <- c(3, 1, 6)
  param <- c(N, odds)
  tvec <- c(0.04, -0.01, -0.02)

  fd <- make_Q_factored(3, rank = 2, seed = 30)

  K3K3_fact <- cgf$.private_api$K3K3operatorABCABC_factored

  got_fact <- K3K3_fact(tvec, param, fd$A, fd$d)
  got_Q    <- cgf$K3K3operatorABCABC(tvec, param, fd$Q)

  expect_equal(as.numeric(got_fact), as.numeric(got_Q), tolerance = 1e-10)
})

# Deferred test:
# Re-enable this once factored methods are implemented broadly enough across
# CGFs that iidReplicatesCGF should consistently delegate in the B == 1 case.
# test_that("iidReplicates delegates factored K3K3 operators to the child when only one block is present", {
#   child <- createCGF(
#     K = function(tvec, parameter_vector) 0 * sum(tvec),
#     K1 = function(tvec, parameter_vector) numeric(length(tvec)) * 0,
#     K2 = function(tvec, parameter_vector) diag(length(tvec)),
#     K3operator = function(tvec, parameter_vector, v1, v2, v3) 0,
#     K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) 0,
#     K3K3operatorAABBCC_factored = function(tvec, parameter_vector, A, d) 123.25,
#     K3K3operatorABCABC_factored = function(tvec, parameter_vector, A, d) 456.5,
#     op_name = "delegatingFactoredToyCGF"
#   )
#
#   wrapped <- iidReplicatesCGF(child, iidReps = "any", block_size = 2)
#   tvec <- c(0.1, -0.2)
#   param <- 1
#   A <- matrix(c(
#     1.0, 0.3,
#     -0.1, 0.9
#   ), nrow = 2, byrow = TRUE)
#   dvec <- c(1.2, 0.7)
#
#   expect_equal(as.numeric(wrapped$.private_api$K3K3operatorAABBCC_factored(tvec, param, A, dvec)), 123.25)
#   expect_equal(as.numeric(wrapped$.private_api$K3K3operatorABCABC_factored(tvec, param, A, dvec)), 456.5)
# })
