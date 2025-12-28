



test_that("randomlyStoppedSumCGF derivatives match numDeriv (nontrivial 3D example)", {
  skip_if_not_installed("numDeriv")

  # Count: N ~ Geometric(tau)  (tau fixed)
  # Summand: X is 3D independent Gamma(shape, rate) (shape/rate fixed)
  count_cgf <- GeometricModelCGF(prob = adaptor(fixed_param = 0.35), iidReps = 1)
  summand_cgf <- GammaModelCGF(shape = adaptor(fixed_param = 2.0),
                               rate  = adaptor(fixed_param = 3.0),
                               iidReps = 3)  # 3-vector of independent Gammas

  cgf <- randomlyStoppedSumCGF(count_cgf = count_cgf,
                               summand_cgf = summand_cgf,
                               iidReps = "any",
                               block_size = 3)

  # No parameters are used (all fixed_param), but keep numeric param for code paths
  param <- c(1.0)

  t0 <- c(0.05, -0.03, 0.01)
  f  <- function(t) cgf$K(t, param)

  # K1 / K2 against numDeriv
  g_fd <- numDeriv::grad(f, t0)
  H_fd <- numDeriv::hessian(f, t0)

  g <- as.numeric(cgf$K1(t0, param))
  H <- as.matrix(cgf$K2(t0, param))

  expect_equal(g, g_fd, tolerance = 5e-6)
  expect_equal(H, H_fd, tolerance = 5e-5)

  # K2operator agrees with x' H y
  set.seed(1)
  x <- rnorm(3)
  y <- rnorm(3)
  ref_xy <- as.numeric(t(x) %*% (H %*% y))
  got_xy <- as.numeric(cgf$K2operator(t0, param, x, y))
  expect_equal(got_xy, ref_xy, tolerance = 1e-8)

  # ---- More stable mixed 3rd/4th checks via numDeriv on g(t) = u1^T K2(t) u2 ----
  set.seed(123)
  u1 <- rnorm(3); u2 <- rnorm(3); u3 <- rnorm(3); u4 <- rnorm(3)

  # scalar function g(t) = u1^T K2(t) u2
  gfun <- function(t) {
    K2t <- as.matrix(cgf$K2(t, param))
    as.numeric(t(u1) %*% (K2t %*% u2))
  }

  # K3(u1,u2,u3) = grad(g)(t0)^T u3
  grad_g <- numDeriv::grad(gfun, t0)
  fd3 <- as.numeric(crossprod(grad_g, u3))
  got3 <- as.numeric(cgf$K3operator(t0, param, u1, u2, u3))
  expect_equal(got3, fd3, tolerance = 2e-4)

  # K4(u1,u2,u3,u4) = u3^T Hess(g)(t0) u4
  Hess_g <- numDeriv::hessian(gfun, t0)
  fd4 <- as.numeric(t(u3) %*% (Hess_g %*% u4))
  got4 <- as.numeric(cgf$K4operator(t0, param, u1, u2, u3, u4))
  expect_equal(got4, fd4, tolerance = 5e-3)
})
