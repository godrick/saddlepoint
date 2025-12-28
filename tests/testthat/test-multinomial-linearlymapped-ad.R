
test_that("compute.spa.negll returns AD gradient/Hessian and matches num-deriv of taped value", {
  skip_on_cran()
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")

  cgfX <- MultinomialCGF

  # Map R^3 -> R^2
  A <- matrix(c(1, 0, 0,
                0, 1, 0), nrow = 2, byrow = TRUE)

  mapped <- linearlyMappedCGF(cgf = cgfX, matrix_A = A, iidReps = "any")

  # Parameter vector: c(N, odds)
  theta0 <- c(12, 2, 3, 5)

  # Observed data for first 2 categories (<= N)
  y <- c(3, 4)

  # AD result (Newton tape)
  res <- compute.spa.negll(
    parameter_vector = theta0,
    observed.data    = y,
    cgf              = mapped,
    gradient         = TRUE,
    hessian          = TRUE,
    tvec_source      = "newton",
    spa_method       = "standard"
  )

  expect_true(is.finite(res$vals))
  expect_true(all(is.finite(res$gradient)))
  expect_true(all(is.finite(res$hessian)))

  # numderiv check against the same taped value function
  # Build a taped evaluator (no derivatives)
  taped_eval <- saddlepoint:::create_spa_taped_fun(
    param_vec     = theta0,
    observed.data = y,
    cgf           = mapped,
    spa_method    = "negll_standard",
    tvec_source   = "newton",
    gradient      = FALSE,
    hessian       = FALSE
  )

  f <- function(th) as.numeric(taped_eval(th)$vals)

  g_fd <- numDeriv::grad(f, theta0)
  H_fd <- numDeriv::hessian(f, theta0)

  expect_equal(as.numeric(res$gradient), as.numeric(g_fd), tolerance = 5e-3)
  expect_equal(as.matrix(res$hessian),  as.matrix(H_fd), tolerance = 5e-2)
})
