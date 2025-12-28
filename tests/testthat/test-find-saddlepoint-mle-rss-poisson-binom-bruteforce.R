
# # Exact log-likelihood for Y = sum_{i=1}^N X_i,
# # where N ~ Poisson(lambda_N) and X_i ~ Binomial(m, p) i.i.d.
# #
# # Conditional on N=n: Y | N=n ~ Binomial(size = n*m, prob = p).
# #
# # P(Y=y) = sum_{n>=ceiling(y/m)} dpois(n;lambda_N) * dbinom(y; n*m, p).
#
# .log_sum_exp <- function(v) {
#   m <- max(v)
#   m + log(sum(exp(v - m)))
# }
#
# .log_pmf_compound_pois_binom <- function(y, p, lambda_N, m, Nmax) {
#   if (!is.finite(p) || p <= 0 || p >= 1) return(-Inf)
#   n_min <- if (y <= 0) 0 else ceiling(y / m)
#   n <- n_min:Nmax
#   # dbinom(y, size=n*m) is defined for size >= y; otherwise it is -Inf.
#   log_terms <- dpois(n, lambda_N, log = TRUE) + dbinom(y, size = n * m, prob = p, log = TRUE)
#   .log_sum_exp(log_terms)
# }
#
# .exact_mle_p_compound_pois_binom <- function(y_vec, lambda_N, m) {
#   y_max <- max(y_vec)
#   # Need n*m >= y_max. ... plus a buffer.
#   Nmax_data <- ceiling(y_max / m) + 40
#   Nmax_tail <- ceiling(lambda_N + 12 * sqrt(lambda_N + 1)) + 40
#   Nmax <- max(Nmax_data, Nmax_tail)
#
#   nll <- function(p) {
#     lp <- vapply(
#       y_vec,
#       function(yy) .log_pmf_compound_pois_binom(yy, p = p, lambda_N = lambda_N, m = m, Nmax = Nmax),
#       numeric(1)
#     )
#     if (any(!is.finite(lp))) return(Inf)
#     -sum(lp)
#   }
#
#   opt <- optimize(nll, interval = c(1e-6, 1 - 1e-6))
#   list(p_hat = opt$minimum, nll = opt$objective, Nmax = Nmax)
# }
#
#
# test_that("RSS Poisson+Binomial(m>1): Exact MLE is close; discrepancy points in the right direction", {
#
#   set.seed(1)
#   B <- 20 # 8
#   lambda_N <- 8 # 1.8
#   m <- 3
#   p_true <- 0.65 # 0.79
#
#   # Simulate from the compound model.
#   N_vec <- rpois(B, lambda = lambda_N)
#   y <- vapply(N_vec, function(n) {
#     if (n == 0) return(0L)
#     rbinom(1, size = n * m, prob = p_true)
#   }, integer(1))
#   # y <- y[y!=0]
#   # B = length(y)
#
#   # Exact MLE for p.
#   exact <- .exact_mle_p_compound_pois_binom(y_vec = y, lambda_N = lambda_N, m = m)
#   p_hat_exact <- exact$p_hat
#
#
#
#
#   count_cgf <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda_N), iidReps = 1)
#   summand_cgf <- BinomialModelCGF(n = adaptor(fixed_param = m), p = adaptor(indices = 1), iidReps = 1)
#   rss_cgf <- randomlyStoppedSumCGF(count_cgf = count_cgf, summand_cgf = summand_cgf, iidReps = B)
#
#   start_theta <- c(0.35)
#   lb_theta <- 0.2
#   ub_theta <- 0.95
#
#   opts <- list(maxeval = 600, xtol_rel = 1e-10, print_level = 0)
#
#   res_spa <- find.saddlepoint.MLE(
#     observed.data = y,
#     cgf = rss_cgf,
#     starting.theta = start_theta,
#     lb.theta = lb_theta,
#     ub.theta = ub_theta,
#     std.error = TRUE,
#     discrepancy = TRUE,
#     opts.user = opts,
#     method = "two_step"
#   )
#
#   p_hat_spa <- res_spa$MLEs.theta
#   disc <- res_spa$discrepancy
#
#   # The SPA MLE should be reasonably close to the brute-force exact MLE.
#   expect_lt(abs(p_hat_spa - p_hat_exact), 0.08)
#
#   # Discrepancy is intended to approximate (exact - spa).
#   delta_exact <- p_hat_exact - p_hat_spa
#   err0 <- abs(delta_exact)
#   err1 <- abs(p_hat_exact - (p_hat_spa + disc))
#
#   # Direction + improvement checks
#   if (abs(disc) > 0) expect_true(delta_exact * disc > 0)
#   expect_true(err1 <= err0 + 1e-12)
#
#   if (err0 > 1e-6) {
#     rel_err <- abs(delta_exact - disc) / err0
#     expect_true(rel_err < 0.5)
#   } else {
#     expect_true(abs(disc) < 1e-6)
#   }
#
#   # Constrained and two_step should agree closely on this 1-parameter model.
#   res_con <- find.saddlepoint.MLE(
#     observed.data = y,
#     cgf = rss_cgf,
#     starting.theta = start_theta,
#     lb.theta = lb_theta,
#     ub.theta = ub_theta,
#     std.error = TRUE,
#     discrepancy = TRUE,
#     opts.user = opts,
#     method = "constrained"
#   )
#   expect_equal(res_con$MLEs.theta, res_spa$MLEs.theta, tolerance = 1e-3)
# })
