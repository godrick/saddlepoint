

test_that("two_step uses nloptr when user inequality constraints are supplied", {

  set.seed(606)
  n <- 10
  y <- rpois(n, lambda = 7.2)

  # Exact unconstrained MLE is mean(y)
  mle_unconstrained <- mean(y)

  # Impose a user constraint: lambda >= 8
  # NLOPT convention in this package: constraints <= 0 is feasible
  user_ineq <- function(theta) {
    c1 <- 8 - theta[1]             # <= 0  <=>  theta[1] >= 8
    J <- matrix(-1, nrow = 1)      # d(c1)/d(theta) = -1
    list(constraints = c1, jacobian = J)
  }

  opts <- list(maxeval = 300, xtol_rel = 1e-10, print_level = 0)

  res <- find.saddlepoint.MLE(
    observed.data = y,
    cgf = PoissonCGF,
    starting.theta = c(1),
    lb.theta = c(1e-8),
    ub.theta = c(Inf),
    method = "two_step",
    user.ineq.constraint.function = user_ineq,
    std.error = TRUE,
    discrepancy = TRUE,
    opts.user = opts
  )

  # two_step should switch to nloptr path
  expect_true(identical(res$method, "two_step"))
  expect_true(identical(res$optimizer, "nloptr"))

  # Constraint should bind when unconstrained MLE < 2
  expect_true(mle_unconstrained < 8)
  expect_true(res$MLEs.theta >= 8 - 1e-6)
  expect_equal(res$MLEs.theta, 8, tolerance = 1e-3)
})
