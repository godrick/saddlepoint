scalarize_operator_value <- function(x) {
  unname(as.numeric(x))
}

testCGF_K3K4ops <- function(cgf,
                            len,
                            parameter_vector,
                            tvec = runif(len),
                            rank = len,
                            B = matrix(rnorm(len * rank), nrow = len, ncol = rank),
                            d = runif(rank)) {
  len <- length(tvec)
  if (!is.matrix(B) || nrow(B) != len || ncol(B) != length(d)) {
    stop("B must be a matrix with length(tvec) rows and length(d) columns.")
  }

  e <- diag(nrow = len)

  naiveK3array <- array(dim = rep(len, times = 3))
  for (i1 in 1:len) {
    for (i2 in 1:len) {
      for (i3 in 1:len) {
        naiveK3array[i1, i2, i3] <- scalarize_operator_value(
          cgf$K3operator(tvec, parameter_vector, e[, i1], e[, i2], e[, i3])
        )
      }
    }
  }

  naiveK4array <- array(dim = rep(len, times = 4))
  for (i1 in 1:len) {
    for (i2 in 1:len) {
      for (i3 in 1:len) {
        for (i4 in 1:len) {
          naiveK4array[i1, i2, i3, i4] <- scalarize_operator_value(
            cgf$K4operator(tvec, parameter_vector, e[, i1], e[, i2], e[, i3], e[, i4])
          )
        }
      }
    }
  }

  Q <- B %*% diag(d, nrow = length(d), ncol = length(d)) %*% t(B)

  naiveK4AABBsummands <- array(dim = rep(len, 4))
  for (i1 in 1:len) {
    for (i2 in 1:len) {
      for (i3 in 1:len) {
        for (i4 in 1:len) {
          naiveK4AABBsummands[i1, i2, i3, i4] <-
            naiveK4array[i1, i2, i3, i4] * Q[i1, i2] * Q[i3, i4]
        }
      }
    }
  }
  naiveK4AABB <- sum(naiveK4AABBsummands)

  naiveK3K3AABBCCsummands <- array(dim = rep(len, 6))
  naiveK3K3ABCABCsummands <- array(dim = rep(len, 6))
  for (i1 in 1:len) {
    for (i2 in 1:len) {
      for (i3 in 1:len) {
        for (j1 in 1:len) {
          for (j2 in 1:len) {
            for (j3 in 1:len) {
              naiveK3K3AABBCCsummands[i1, i2, i3, j1, j2, j3] <-
                naiveK3array[i1, i2, i3] * naiveK3array[j1, j2, j3] *
                Q[i1, i2] * Q[i3, j1] * Q[j2, j3]
              naiveK3K3ABCABCsummands[i1, i2, i3, j1, j2, j3] <-
                naiveK3array[i1, i2, i3] * naiveK3array[j1, j2, j3] *
                Q[i1, j1] * Q[i2, j2] * Q[i3, j3]
            }
          }
        }
      }
    }
  }
  naiveK3K3AABBCC <- sum(naiveK3K3AABBCCsummands)
  naiveK3K3ABCABC <- sum(naiveK3K3ABCABCsummands)

  compare_values <- function(actual, expected) {
    all.equal(scalarize_operator_value(actual), scalarize_operator_value(expected))
  }

  res <- list(
    all_passed = NULL,
    K4agrees = compare_values(cgf$K4operatorAABB(tvec, parameter_vector, Q), naiveK4AABB),
    K3K3AABBCCagrees = compare_values(cgf$K3K3operatorAABBCC(tvec, parameter_vector, Q), naiveK3K3AABBCC),
    K3K3ABCABCagrees = compare_values(cgf$K3K3operatorABCABC(tvec, parameter_vector, Q), naiveK3K3ABCABC),
    K4factored_agrees = compare_values(cgf$.private_api$K4operatorAABB_factored(tvec, parameter_vector, B, d), naiveK4AABB),
    K3K3AABBCCfactored_agrees = compare_values(cgf$.private_api$K3K3operatorAABBCC_factored(tvec, parameter_vector, B, d), naiveK3K3AABBCC),
    K3K3ABCABCfactored_agrees = compare_values(cgf$.private_api$K3K3operatorABCABC_factored(tvec, parameter_vector, B, d), naiveK3K3ABCABC),
    K4_consistent = compare_values(
      cgf$K4operatorAABB(tvec, parameter_vector, Q),
      cgf$.private_api$K4operatorAABB_factored(tvec, parameter_vector, B, d)
    ),
    K3K3AABBCC_consistent = compare_values(
      cgf$K3K3operatorAABBCC(tvec, parameter_vector, Q),
      cgf$.private_api$K3K3operatorAABBCC_factored(tvec, parameter_vector, B, d)
    ),
    K3K3ABCABC_consistent = compare_values(
      cgf$K3K3operatorABCABC(tvec, parameter_vector, Q),
      cgf$.private_api$K3K3operatorABCABC_factored(tvec, parameter_vector, B, d)
    )
  )
  res[[1]] <- all(vapply(res[-1], isTRUE, logical(1)))
  attributes(res)$cgf <- cgf
  attributes(res)$call_history <- cgf$call_history
  res
}

summarize_testCGF_K3K4ops <- function(res) {
  checks <- names(res)[-1]
  data.frame(
    check = checks,
    result = vapply(res[-1], function(x) if (isTRUE(x)) "PASS" else as.character(x), character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_testCGF_K3K4ops_cases <- function() {
  N <- 4
  M <- 3
  theta_pois <- 1.5
  theta_gamma <- c(shape = 3, rate = 1.5)
  A_map <- matrix(c(
     0.8, -0.2,  0.4,  0.1,
    -0.3,  1.0,  0.2, -0.4,
     0.1,  0.3,  0.9,  0.2
  ), nrow = M, byrow = TRUE)

  list(
    poisson_base = list(
      cgf = PoissonCGF,
      len = N,
      parameter_vector = theta_pois,
      tvec = c(0.10, 0.20, 0.05, 0.15),
      rank = N,
      B = matrix(c(
        1.0,  0.1, -0.1,  0.2,
        0.0,  0.9,  0.3, -0.2,
        0.2, -0.1,  1.1,  0.0,
       -0.1,  0.2,  0.0,  0.8
      ), nrow = N, byrow = TRUE),
      d = c(1.2, 0.9, 1.1, 0.8)
    ),
    concat_poisson = list(
      cgf = concatenationCGF(rep(list(PoissonCGF), times = N)),
      len = N,
      parameter_vector = theta_pois,
      tvec = c(0.10, 0.20, 0.05, 0.15),
      rank = N,
      B = matrix(c(
        1.0,  0.1, -0.1,  0.2,
        0.0,  0.9,  0.3, -0.2,
        0.2, -0.1,  1.1,  0.0,
       -0.1,  0.2,  0.0,  0.8
      ), nrow = N, byrow = TRUE),
      d = c(1.2, 0.9, 1.1, 0.8)
    ),
    linmap_poisson = list(
      cgf = linearlyMappedCGF(PoissonCGF, A_map),
      len = M,
      parameter_vector = theta_pois,
      tvec = c(0.08, 0.16, 0.05),
      rank = M,
      B = matrix(c(
        0.9, -0.2,  0.3,
        0.1,  1.0, -0.1,
       -0.2,  0.1,  0.8
      ), nrow = M, byrow = TRUE),
      d = c(1.1, 0.7, 0.9)
    ),
    linmap_concat_poisson = list(
      cgf = linearlyMappedCGF(concatenationCGF(rep(list(PoissonCGF), times = N)), A_map),
      len = M,
      parameter_vector = theta_pois,
      tvec = c(0.08, 0.16, 0.05),
      rank = M,
      B = matrix(c(
        0.9, -0.2,  0.3,
        0.1,  1.0, -0.1,
       -0.2,  0.1,  0.8
      ), nrow = M, byrow = TRUE),
      d = c(1.1, 0.7, 0.9)
    ),
    gamma_base = list(
      cgf = GammaCGF,
      len = N,
      parameter_vector = theta_gamma,
      tvec = c(0.10, 0.20, 0.05, 0.15),
      rank = N,
      B = matrix(c(
        1.0,  0.1, -0.1,  0.2,
        0.0,  0.9,  0.3, -0.2,
        0.2, -0.1,  1.1,  0.0,
       -0.1,  0.2,  0.0,  0.8
      ), nrow = N, byrow = TRUE),
      d = c(1.2, 0.9, 1.1, 0.8)
    ),
    linmap_gamma = list(
      cgf = linearlyMappedCGF(GammaCGF, A_map),
      len = M,
      parameter_vector = theta_gamma,
      tvec = c(0.08, 0.16, 0.05),
      rank = M,
      B = matrix(c(
        0.9, -0.2,  0.3,
        0.1,  1.0, -0.1,
       -0.2,  0.1,  0.8
      ), nrow = M, byrow = TRUE),
      d = c(1.1, 0.7, 0.9)
    ),
    iidrep_scalar_poisson = list(
      cgf = iidReplicatesCGF(PoissonCGF, iidReps = 2, block_size = 1),
      len = 2,
      parameter_vector = theta_pois,
      tvec = c(0.10, 0.20),
      rank = 2,
      B = matrix(c(
        1.0,  0.4,
       -0.2,  0.9
      ), nrow = 2, byrow = TRUE),
      d = c(1.0, 0.7)
    ),
    iidrep_linmap_gamma = list(
      cgf = iidReplicatesCGF(linearlyMappedCGF(GammaCGF, A_map), iidReps = 2, block_size = M),
      len = 2 * M,
      parameter_vector = theta_gamma,
      tvec = c(0.04, 0.08, 0.03, 0.05, 0.07, 0.02),
      rank = 4,
      B = matrix(c(
        1.0,  0.2, -0.1,  0.0,
        0.1,  0.9,  0.3, -0.2,
       -0.2,  0.1,  0.8,  0.4,
        0.0, -0.3,  0.2,  1.0,
        0.2,  0.0,  0.5, -0.1,
       -0.1,  0.4,  0.0,  0.9
      ), nrow = 2 * M, byrow = TRUE),
      d = c(1.0, 0.8, 0.9, 0.7)
    )
  )
}

expect_testCGF_K3K4ops_passes <- function(case_name, res) {
  failed_checks <- names(res)[-1][!vapply(res[-1], isTRUE, logical(1))]
  expect(
    res$all_passed,
    paste0(
      "testCGF_K3K4ops failed for ",
      case_name,
      ": ",
      paste(failed_checks, collapse = ", ")
    )
  )
}
