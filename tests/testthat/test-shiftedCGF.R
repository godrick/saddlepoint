
test_that("shiftedCGF: K/K1/K2 identities for a deterministic shift", {
  base <- PoissonCGF
  b <- 3
  cgf <- shiftedCGF(base, shift = b)

  lambda <- 2
  tvec <- c(0.10, -0.20, 0.00, 0.30)

  expect_equal(
    cgf$K(tvec, lambda),
    base$K(tvec, lambda) + sum(tvec * b),
    tolerance = 1e-12
  )
  expect_equal(
    cgf$K1(tvec, lambda),
    base$K1(tvec, lambda) + rep.int(b, length(tvec)),
    tolerance = 1e-12
  )
  expect_equal(
    as.matrix(cgf$K2(tvec, lambda)),
    as.matrix(base$K2(tvec, lambda)),
    tolerance = 1e-12
  )

  te0 <- base$.get_private_method("tilting_exponent")(tvec, lambda)
  te1 <- cgf$.get_private_method("tilting_exponent")(tvec, lambda)
  expect_equal(as.numeric(te1), as.numeric(te0), tolerance = 1e-12)
})


test_that("shiftedCGF: analytic saddlepoint t-hat shifts y -> y - b", {
  base <- NormalCGF
  b <- 1.5
  cgf <- shiftedCGF(base, shift = b)

  theta <- c(mu = 0.2, sigma = 1.1)
  y <- c(0.5, -0.1, 1.2)

  t_hat0 <- base$analytic_tvec_hat(y - b, theta)
  t_hat1 <- cgf$analytic_tvec_hat(y, theta)

  expect_equal(t_hat1, t_hat0, tolerance = 1e-12)
  expect_equal(cgf$K1(t_hat1, theta), y, tolerance = 1e-12)
})


test_that("shiftedCGF: rsim shifts draws when the base CGF can simulate", {
  base <- PoissonCGF
  b <- 3
  cgf <- shiftedCGF(base, shift = b)

  lambda <- 2
  iidReps <- 10L

  set.seed(123)
  x <- base$rsim(iidReps = iidReps, parameter_vector = lambda, drop = FALSE)
  set.seed(123)
  y <- cgf$rsim(iidReps = iidReps, parameter_vector = lambda, drop = FALSE)

  expect_equal(y, x + b)
})


test_that("shiftedCGF: rsim supports vector shifts", {
  base <- PoissonCGF
  b <- c(1, 2)
  cgf <- shiftedCGF(base, shift = b)

  lambda <- c(2, 5)
  iidReps <- 7L

  set.seed(1)
  x <- base$rsim(iidReps = iidReps, parameter_vector = lambda, drop = FALSE)
  set.seed(1)
  y <- cgf$rsim(iidReps = iidReps, parameter_vector = lambda, drop = FALSE)

  expect_equal(y, x + b)
})

