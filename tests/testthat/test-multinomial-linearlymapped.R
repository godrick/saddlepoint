
test_that("Mapping identities hold for a single block", {
  cgfX <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)

  # Map R^3 -> R^2 by dropping the 3rd category
  A <- matrix(c(1, 0, 0,
                0, 1, 0), nrow = 2, byrow = TRUE)

  mapped <- linearlyMappedCGF(cgf = cgfX, matrix_A = A, iidReps = "any")

  tY <- c(0.10, -0.05)  # length 2
  tX <- as.vector(t(A) %*% tY)  # length 3

  # K
  expect_equal(
    as.numeric(mapped$K(tY, param)),
    as.numeric(cgfX$K(tX, param)),
    tolerance = 1e-12
  )

  # K1
  lhs_K1 <- as.numeric(mapped$K1(tY, param))
  rhs_K1 <- as.numeric(A %*% cgfX$K1(tX, param))
  expect_equal(lhs_K1, rhs_K1, tolerance = 1e-12)

  # K2
  lhs_K2 <- as.matrix(mapped$K2(tY, param))
  rhs_K2 <- as.matrix(A %*% cgfX$K2(tX, param) %*% t(A))
  expect_equal(lhs_K2, rhs_K2, tolerance = 1e-12)

  # K3operator
  set.seed(1)
  u1 <- rnorm(2); u2 <- rnorm(2); u3 <- rnorm(2)
  lhs_K3 <- mapped$K3operator(tY, param, u1, u2, u3)
  rhs_K3 <- cgfX$K3operator(tX, param, as.vector(t(A) %*% u1), as.vector(t(A) %*% u2), as.vector(t(A) %*% u3))
  expect_equal(as.numeric(lhs_K3), as.numeric(rhs_K3), tolerance = 1e-10)

  # K4operator
  set.seed(2)
  v1 <- rnorm(2); v2 <- rnorm(2); v3 <- rnorm(2); v4 <- rnorm(2)
  lhs_K4 <- mapped$K4operator(tY, param, v1, v2, v3, v4)
  rhs_K4 <- cgfX$K4operator(tX, param,
                           as.vector(t(A) %*% v1), as.vector(t(A) %*% v2),
                           as.vector(t(A) %*% v3), as.vector(t(A) %*% v4))
  expect_equal(as.numeric(lhs_K4), as.numeric(rhs_K4), tolerance = 1e-10)
})
