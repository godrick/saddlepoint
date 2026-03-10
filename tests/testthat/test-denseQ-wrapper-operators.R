make_toy_diag_cgf <- function(k2, k3, k4) {
  d <- length(k2)
  stopifnot(length(k3) == d, length(k4) == d)

  createCGF(
    K = function(tvec, parameter_vector) {
      sum(0.5 * k2 * tvec^2 + (k3 / 6) * tvec^3 + (k4 / 24) * tvec^4)
    },
    K1 = function(tvec, parameter_vector) {
      k2 * tvec + 0.5 * k3 * tvec^2 + (k4 / 6) * tvec^3
    },
    K2 = function(tvec, parameter_vector) {
      diag(k2 + k3 * tvec + 0.5 * k4 * tvec^2, nrow = d)
    },
    K3operator = function(tvec, parameter_vector, v1, v2, v3) {
      coeff <- k3 + k4 * tvec
      sum(coeff * v1 * v2 * v3)
    },
    K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) {
      sum(k4 * v1 * v2 * v3 * v4)
    },
    op_name = "toyDiagCGF"
  )
}

basis_vector <- function(n, i) {
  v <- numeric(n)
  v[i] <- 1
  v
}

naive_K3_tensor <- function(cgf, tvec, param) {
  d <- length(tvec)
  out <- array(0, dim = c(d, d, d))

  for (i in seq_len(d)) {
    ei <- basis_vector(d, i)
    for (j in seq_len(d)) {
      ej <- basis_vector(d, j)
      for (k in seq_len(d)) {
        ek <- basis_vector(d, k)
        out[i, j, k] <- as.numeric(cgf$K3operator(tvec, param, ei, ej, ek))
      }
    }
  }

  out
}

naive_K4_tensor <- function(cgf, tvec, param) {
  d <- length(tvec)
  out <- array(0, dim = c(d, d, d, d))

  for (i in seq_len(d)) {
    ei <- basis_vector(d, i)
    for (j in seq_len(d)) {
      ej <- basis_vector(d, j)
      for (k in seq_len(d)) {
        ek <- basis_vector(d, k)
        for (l in seq_len(d)) {
          el <- basis_vector(d, l)
          out[i, j, k, l] <- as.numeric(cgf$K4operator(tvec, param, ei, ej, ek, el))
        }
      }
    }
  }

  out
}

naive_K4operatorAABB <- function(K4_tensor, Q) {
  d <- dim(K4_tensor)[1]
  total <- 0

  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      for (k in seq_len(d)) {
        for (l in seq_len(d)) {
          total <- total + K4_tensor[i, j, k, l] * Q[i, j] * Q[k, l]
        }
      }
    }
  }

  total
}

naive_K3K3operatorAABBCC <- function(K3_tensor, Q) {
  d <- dim(K3_tensor)[1]
  total <- 0

  for (a in seq_len(d)) {
    for (b in seq_len(d)) {
      for (c in seq_len(d)) {
        left_val <- K3_tensor[a, b, c]
        for (i in seq_len(d)) {
          for (j in seq_len(d)) {
            for (k in seq_len(d)) {
              total <- total + left_val * K3_tensor[i, j, k] * Q[a, b] * Q[i, j] * Q[c, k]
            }
          }
        }
      }
    }
  }

  total
}

naive_K3K3operatorABCABC <- function(K3_tensor, Q) {
  d <- dim(K3_tensor)[1]
  total <- 0

  for (a in seq_len(d)) {
    for (b in seq_len(d)) {
      for (c in seq_len(d)) {
        left_val <- K3_tensor[a, b, c]
        for (i in seq_len(d)) {
          for (j in seq_len(d)) {
            for (k in seq_len(d)) {
              total <- total + left_val * K3_tensor[i, j, k] * Q[a, i] * Q[b, j] * Q[c, k]
            }
          }
        }
      }
    }
  }

  total
}

naive_operator_bundle <- function(cgf, tvec, param, Q) {
  K3_tensor <- naive_K3_tensor(cgf, tvec, param)
  K4_tensor <- naive_K4_tensor(cgf, tvec, param)

  list(
    K4operatorAABB = naive_K4operatorAABB(K4_tensor, Q),
    K3K3operatorAABBCC = naive_K3K3operatorAABBCC(K3_tensor, Q),
    K3K3operatorABCABC = naive_K3K3operatorABCABC(K3_tensor, Q)
  )
}

expect_public_bundle_matches_naive <- function(cgf, tvec, param, Q, tol = 1e-10) {
  ref <- naive_operator_bundle(cgf, tvec, param, Q)

  expect_equal(
    as.numeric(cgf$K4operatorAABB(tvec, param, Q)),
    as.numeric(ref$K4operatorAABB),
    tolerance = tol
  )
  expect_equal(
    as.numeric(cgf$K3K3operatorAABBCC(tvec, param, Q)),
    as.numeric(ref$K3K3operatorAABBCC),
    tolerance = tol
  )
  expect_equal(
    as.numeric(cgf$K3K3operatorABCABC(tvec, param, Q)),
    as.numeric(ref$K3K3operatorABCABC),
    tolerance = tol
  )
}

expect_factored_bundle_matches_naive <- function(cgf, tvec, param, B, dvec, tol = 1e-10) {
  Q <- B %*% diag(dvec, nrow = length(dvec)) %*% t(B)
  ref <- naive_operator_bundle(cgf, tvec, param, Q)

  expect_equal(
    as.numeric(cgf$.private_api$K4operatorAABB_factored(tvec, param, B, dvec)),
    as.numeric(ref$K4operatorAABB),
    tolerance = tol
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorAABBCC_factored(tvec, param, B, dvec)),
    as.numeric(ref$K3K3operatorAABBCC),
    tolerance = tol
  )
  expect_equal(
    as.numeric(cgf$.private_api$K3K3operatorABCABC_factored(tvec, param, B, dvec)),
    as.numeric(ref$K3K3operatorABCABC),
    tolerance = tol
  )
}

test_that("iidReplicates dense-Q and factored K3K3 operators match naive contractions", {
  toy2 <- make_toy_diag_cgf(
    k2 = c(1.4, 1.1),
    k3 = c(0.7, -0.5),
    k4 = c(0.6, 0.9)
  )
  cgf <- iidReplicatesCGF(toy2, iidReps = 2, block_size = 2)

  tvec <- c(0.10, -0.15, -0.05, 0.20)
  param <- 1
  R <- matrix(c(
    1.0,  0.1, -0.2,  0.0,
    0.2,  0.9,  0.1, -0.1,
   -0.1,  0.3,  0.8,  0.2,
    0.0, -0.2,  0.4,  1.1
  ), nrow = 4, byrow = TRUE)
  Q <- crossprod(R)

  expect_public_bundle_matches_naive(cgf, tvec, param, Q)

  B <- matrix(c(
    1.0,  0.1, -0.2,
    0.2,  0.9,  0.1,
   -0.1,  0.3,  0.8,
    0.0, -0.2,  0.4
  ), nrow = 4, byrow = TRUE)
  dvec <- c(1.3, 0.8, 1.1)

  expect_factored_bundle_matches_naive(cgf, tvec, param, B, dvec)
})

test_that("concatenation dense-Q and factored K3K3 operators match naive contractions", {
  toy2 <- make_toy_diag_cgf(
    k2 = c(1.4, 1.1),
    k3 = c(0.7, -0.5),
    k4 = c(0.6, 0.9)
  )
  toy1 <- make_toy_diag_cgf(
    k2 = c(1.2, 0.8),
    k3 = c(0.9, 0.4),
    k4 = c(0.4, 0.7)
  )
  cgf <- .concatenationCGF_internal(
    cgf_list = list(toy2, toy1),
    component_dims = c(2L, 2L)
  )

  tvec <- c(0.12, -0.08, 0.07, -0.04)
  param <- 1
  R <- matrix(c(
    1.0,  0.1, -0.1,  0.0,
    0.0,  0.9,  0.2, -0.1,
    0.2, -0.1,  1.0,  0.3,
   -0.1,  0.2,  0.0,  0.8
  ), nrow = 4, byrow = TRUE)
  Q <- crossprod(R)

  expect_public_bundle_matches_naive(cgf, tvec, param, Q)

  B <- matrix(c(
    1.0,  0.1, -0.2,
    0.2,  0.9,  0.1,
   -0.1,  0.4,  0.8,
    0.0, -0.2,  0.3
  ), nrow = 4, byrow = TRUE)
  dvec <- c(1.1, 0.7, 0.9)

  expect_factored_bundle_matches_naive(cgf, tvec, param, B, dvec)
})

test_that("linearlyMapped over iidReplicates inherits the corrected dense-Q operators", {
  toy2 <- make_toy_diag_cgf(
    k2 = c(1.4, 1.1),
    k3 = c(0.7, -0.5),
    k4 = c(0.6, 0.9)
  )
  child <- iidReplicatesCGF(toy2, iidReps = 2, block_size = 2)
  A_map <- matrix(c(
    1.0,  0.0,  0.2, -0.1,
   -0.3,  1.0,  0.4,  0.0,
    0.1,  0.2, -0.2,  1.1
  ), nrow = 3, byrow = TRUE)
  cgf <- linearlyMappedCGF(cgf = child, matrix_A = A_map, iidReps = 1)

  tvec <- c(0.06, -0.03, 0.09)
  param <- 1
  R <- matrix(c(
    1.0,  0.2, -0.1,
    0.1,  1.0,  0.0,
   -0.1,  0.1,  0.9
  ), nrow = 3, byrow = TRUE)
  Q <- crossprod(R)

  expect_public_bundle_matches_naive(cgf, tvec, param, Q)
})

test_that("iidReplicates delegates dense K3K3 operators to the child when only one block is present", {
  child <- createCGF(
    K = function(tvec, parameter_vector) sum(tvec) * 0,
    K1 = function(tvec, parameter_vector) numeric(length(tvec)) * 0,
    K2 = function(tvec, parameter_vector) diag(length(tvec)),
    K3operator = function(tvec, parameter_vector, v1, v2, v3) 0,
    K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) 0,
    K3K3operatorAABBCC = function(tvec, parameter_vector, Q) 123.25,
    K3K3operatorABCABC = function(tvec, parameter_vector, Q) 456.5,
    op_name = "delegatingToyCGF"
  )

  wrapped <- iidReplicatesCGF(child, iidReps = "any", block_size = 2)
  tvec <- c(0.1, -0.2)
  param <- 1
  Q <- matrix(c(1.0, 0.3,
                0.3, 0.9), nrow = 2, byrow = TRUE)

  expect_equal(as.numeric(wrapped$K3K3operatorAABBCC(tvec, param, Q)), 123.25)
  expect_equal(as.numeric(wrapped$K3K3operatorABCABC(tvec, param, Q)), 456.5)
})
