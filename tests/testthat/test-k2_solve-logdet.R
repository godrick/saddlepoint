
test_that("If implemented, K2_solve and logdetK2 agree with solve/determinant", {
  cgf <- MultinomialCGF

  # Skip unless these hooks exist
  if (is.null(cgf$K2_solve) || is.null(cgf$logdetK2)) {
    skip("K2_solve/logdetK2 not implemented on CGF objects yet")
  }

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
