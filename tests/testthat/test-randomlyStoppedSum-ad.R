
test_that("compute.spa.negll returns AD gradient/Hessian and matches numDeriv on taped value", {
  skip_on_cran()
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")

  # A non-degenerate (identifiable) 2-parameter scalar RSS:
  # N ~ Poisson(lambda), X ~ Binomial(n=2, p)
  # K_Y(t) = lambda * ( (1-p+p e^t)^2 - 1 )
  count <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summ  <- BinomialModelCGF(n = adaptor(fixed_param = 2), p = adaptor(indices = 2), iidReps = 1)

  B <- 5L
  rss <- randomlyStoppedSumCGF(count, summ, iidReps = B) # block_size inferred as 1

  theta0 <- c(lambda = 1.1, p = 0.35)

  # Some plausible observed counts (scalar per replicate)
  y <- c(1, 1, 2, 1, 3)

  res <- compute.spa.negll(
    parameter_vector = theta0,
    observed.data    = y,
    cgf              = rss,
    gradient         = TRUE,
    hessian          = TRUE,
    tvec_source      = "newton",
    spa_method       = "standard"
  )

  expect_true(is.finite(res$vals))
  expect_true(all(is.finite(res$gradient)))
  expect_true(all(is.finite(res$hessian)))

  taped_eval <- saddlepoint:::create_spa_taped_fun(
    param_vec     = theta0,
    observed.data = y,
    cgf           = rss,
    spa_method    = "negll_standard",
    tvec_source   = "newton",
    gradient      = FALSE,
    hessian       = FALSE
  )

  f <- function(th) as.numeric(taped_eval(th)$vals)

  g_fd <- numDeriv::grad(func = f, x = theta0)
  H_fd <- numDeriv::hessian(func = f, x = theta0)

  expect_equal(as.numeric(res$gradient), as.numeric(g_fd), tolerance = 5e-3)
  expect_equal(as.matrix(res$hessian),  as.matrix(H_fd), tolerance = 5e-2)
})
