#
# fd_mixed3 <- function(f, t0, u1, u2, u3, h = 1e-5) {
#   # Mixed third derivative along directions u1,u2,u3 using 8-point central difference
#   stopifnot(length(t0) == length(u1), length(u1) == length(u2), length(u2) == length(u3))
#   signs <- c(-1, 1)
#   acc <- 0
#   for (s1 in signs) for (s2 in signs) for (s3 in signs) {
#     acc <- acc + (s1 * s2 * s3) * f(t0 + h * (s1 * u1 + s2 * u2 + s3 * u3))
#   }
#   acc / (8 * h^3)
# }
#
# fd_mixed4 <- function(f, t0, u1, u2, u3, u4, h = 1e-4) {
#   # Mixed fourth derivative along directions u1,u2,u3,u4 using 16-point central difference
#   stopifnot(length(t0) == length(u1), length(u1) == length(u2), length(u2) == length(u3), length(u3) == length(u4))
#   signs <- c(-1, 1)
#   acc <- 0
#   for (s1 in signs) for (s2 in signs) for (s3 in signs) for (s4 in signs) {
#     acc <- acc + (s1 * s2 * s3 * s4) * f(t0 + h * (s1 * u1 + s2 * u2 + s3 * u3 + s4 * u4))
#   }
#   acc / (16 * h^4)
# }
#
#
# test_that("K/K1/K2 match closed form", {
#   cgf <- MultinomialCGF
#
#   N <- 10
#   odds <- c(2, 3, 5)
#   p <- odds / sum(odds)
#   param <- c(N, odds)
#
#   tvec <- c(0.10, -0.20, 0.05)
#
#   # closed form
#   denom <- sum(p * exp(tvec))
#   K_ref <- N * log(denom)
#   v <- (p * exp(tvec)) / denom
#   K1_ref <- N * v
#   K2_ref <- N * (diag(v) - outer(v, v))
#
#   expect_equal(cgf$K(tvec, param), K_ref, tolerance = 1e-12)
#   expect_equal(as.numeric(cgf$K1(tvec, param)), as.numeric(K1_ref), tolerance = 1e-12)
#   expect_equal(as.matrix(cgf$K2(tvec, param)), as.matrix(K2_ref), tolerance = 1e-12)
#
#   # sanity
#   expect_equal(t(cgf$K2(tvec, param)), cgf$K2(tvec, param), tolerance = 1e-12)
# })
#
#
# test_that("K2operator and K2operatorAK2AT agree with K2", {
#   cgf <- MultinomialCGF
#
#   N <- 12
#   odds <- c(0.2, 0.3, 0.5)  # probabilities (also valid)
#   param <- c(N, odds)
#
#   tvec <- c(0.12, -0.07, 0.02)
#   K2 <- as.matrix(cgf$K2(tvec, param))
#
#   set.seed(1)
#   x <- rnorm(3)
#   y <- rnorm(3)
#
#   ref_xy <- as.numeric(t(x) %*% (K2 %*% y))
#   got_xy <- as.numeric(cgf$K2operator(tvec, param, x, y))
#   expect_equal(got_xy, ref_xy, tolerance = 1e-10)
#
#   A <- matrix(rnorm(2 * 3), nrow = 2)
#   ref_A <- A %*% K2 %*% t(A)
#   got_A <- cgf$K2operatorAK2AT(tvec, param, A)
#   expect_equal(as.matrix(got_A), as.matrix(ref_A), tolerance = 1e-10)
# })
#
#
# test_that("K3operator matches finite-difference mixed third derivative", {
#   cgf <- MultinomialCGF
#
#   N <- 10
#   odds <- c(2, 3, 5)
#   param <- c(N, odds)
#
#   tvec <- c(0.05, -0.08, 0.02)
#
#   set.seed(2)
#   u1 <- rnorm(3)
#   u2 <- rnorm(3)
#   u3 <- rnorm(3)
#
#   f <- function(t) cgf$K(t, param)
#
#   fd <- fd_mixed3(f, tvec, u1, u2, u3, h = 1e-5)
#   got <- cgf$K3operator(tvec, param, u1, u2, u3)
#
#   expect_equal(as.numeric(got), as.numeric(fd), tolerance = 1e-5)
#
#   # symmetry sanity: swap u2 and u3
#   got2 <- cgf$K3operator(tvec, param, u1, u3, u2)
#   expect_equal(as.numeric(got2), as.numeric(got), tolerance = 1e-12)
# })
#
#
# test_that("K4operator matches finite-difference mixed fourth derivative", {
#   cgf <- MultinomialCGF
#
#   N <- 10
#   odds <- c(2, 3, 5)
#   param <- c(N, odds)
#
#   tvec <- c(0.03, -0.04, 0.01)
#
#   set.seed(3)
#   u1 <- rnorm(3)
#   u2 <- rnorm(3)
#   u3 <- rnorm(3)
#   u4 <- rnorm(3)
#
#   f <- function(t) cgf$K(t, param)
#
#   fd <- fd_mixed4(f, tvec, u1, u2, u3, u4, h = 5e-4)
#   got <- cgf$K4operator(tvec, param, u1, u2, u3, u4)
#
#   expect_equal(as.numeric(got), as.numeric(fd), tolerance = 1e-4)
#
#   # symmetry sanity: swap u1 and u2
#   got2 <- cgf$K4operator(tvec, param, u2, u1, u3, u4)
#   expect_equal(as.numeric(got2), as.numeric(got), tolerance = 1e-12)
# })




test_that("K/K1/K2 match closed form", {
  cgf <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  p <- odds / sum(odds)
  param <- c(N, odds)

  tvec <- c(0.10, -0.20, 0.05)

  # closed form
  denom <- sum(p * exp(tvec))
  K_ref <- N * log(denom)
  v <- (p * exp(tvec)) / denom
  K1_ref <- N * v
  K2_ref <- N * (diag(v) - outer(v, v))

  expect_equal(cgf$K(tvec, param), K_ref, tolerance = 1e-12)
  expect_equal(as.numeric(cgf$K1(tvec, param)), as.numeric(K1_ref), tolerance = 1e-12)
  expect_equal(as.matrix(cgf$K2(tvec, param)), as.matrix(K2_ref), tolerance = 1e-12)

  # sanity: symmetry
  expect_equal(t(cgf$K2(tvec, param)), cgf$K2(tvec, param), tolerance = 1e-12)
})


test_that("K2operator and K2operatorAK2AT agree with K2", {
  cgf <- MultinomialCGF

  N <- 12
  odds <- c(0.2, 0.3, 0.5)  # probabilities (also valid)
  param <- c(N, odds)

  tvec <- c(0.12, -0.07, 0.02)
  K2 <- as.matrix(cgf$K2(tvec, param))

  set.seed(1)
  x <- rnorm(3)
  y <- rnorm(3)

  ref_xy <- as.numeric(t(x) %*% (K2 %*% y))
  got_xy <- as.numeric(cgf$K2operator(tvec, param, x, y))
  expect_equal(got_xy, ref_xy, tolerance = 1e-10)

  # A is r x d; K2operatorAK2AT should return r x r
  A <- matrix(rnorm(2 * 3), nrow = 2)
  ref_A <- A %*% K2 %*% t(A)
  got_A <- cgf$K2operatorAK2AT(tvec, param, A)
  expect_equal(as.matrix(got_A), as.matrix(ref_A), tolerance = 1e-10)
})


test_that("K3operator matches numeric derivative via numDeriv::grad", {
  testthat::skip_if_not_installed("numDeriv")

  cgf <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)

  tvec <- c(0.05, -0.08, 0.02)

  set.seed(2)
  u1 <- rnorm(3)
  u2 <- rnorm(3)
  u3 <- rnorm(3)

  # g(t) = u1^T K2(t) u2
  g <- function(t) {
    K2 <- as.matrix(cgf$K2(t, param))
    as.numeric(t(u1) %*% (K2 %*% u2))
  }

  grad_g <- numDeriv::grad(g, tvec)
  fd <- as.numeric(crossprod(grad_g, u3))  # directional derivative along u3

  got <- cgf$K3operator(tvec, param, u1, u2, u3)

  expect_equal(as.numeric(got), as.numeric(fd), tolerance = 1e-5)

  # symmetry sanity: swap u2 and u3
  got2 <- cgf$K3operator(tvec, param, u1, u3, u2)
  expect_equal(as.numeric(got2), as.numeric(got), tolerance = 1e-12)
})


test_that("K4operator matches numeric derivative via numDeriv::hessian", {
  testthat::skip_if_not_installed("numDeriv")

  cgf <- MultinomialCGF

  N <- 10
  odds <- c(2, 3, 5)
  param <- c(N, odds)

  tvec <- c(0.03, -0.04, 0.01)

  set.seed(3)
  u1 <- rnorm(3)
  u2 <- rnorm(3)
  u3 <- rnorm(3)
  u4 <- rnorm(3)

  # g(t) = u1^T K2(t) u2
  g <- function(t) {
    K2 <- as.matrix(cgf$K2(t, param))
    as.numeric(t(u1) %*% (K2 %*% u2))
  }

  H <- numDeriv::hessian(g, tvec)
  fd <- as.numeric(t(u3) %*% (H %*% u4))  # bilinear form u3^T Hess(g) u4

  got <- cgf$K4operator(tvec, param, u1, u2, u3, u4)

  expect_equal(as.numeric(got), as.numeric(fd), tolerance = 1e-4)

  # symmetry sanity: swap u1 and u2
  got2 <- cgf$K4operator(tvec, param, u2, u1, u3, u4)
  expect_equal(as.numeric(got2), as.numeric(got), tolerance = 1e-12)
})
