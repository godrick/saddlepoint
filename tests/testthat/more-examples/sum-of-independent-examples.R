# ## ============================================================
# ## Compare old vs new sumOfIndependentCGF + examples
# ## ============================================================
#
# ## ------------------------------------------------------------
# ## A "slow" baseline that does NOT override factored operators
# ## ------------------------------------------------------------
# .sumOfIndependentCGF_internal_slow <- function(cgf_list, ...) {
#
#   if (!is.list(cgf_list) || length(cgf_list) == 0L) {
#     stop("'cgf_list' must be a non-empty list of CGF objects.")
#   }
#   if (any(vapply(cgf_list, function(x) !inherits(x, "CGF"), logical(1)))) {
#     stop("Every element of 'cgf_list' must be of class 'CGF'.")
#   }
#
#   K_list    <- lapply(cgf_list, function(cg) cg$K)
#   K1_list   <- lapply(cgf_list, function(cg) cg$K1)
#   K2_list   <- lapply(cgf_list, function(cg) cg$K2)
#   K3_list   <- lapply(cgf_list, function(cg) cg$K3operator)
#   K4_list   <- lapply(cgf_list, function(cg) cg$K4operator)
#
#   tilting_list <- lapply(cgf_list, function(cg) cg$.get_private_method("tilting_exponent"))
#   ineq_list    <- lapply(cgf_list, function(cg) cg$ineq_constraint)
#
#   Kfun <- function(tvec, param) {
#     total <- 0 * param[1]
#     for (f in K_list) total <- total + f(tvec, param)
#     total
#   }
#   K1fun <- function(tvec, param) {
#     out <- numeric(length(tvec)) * param[1]
#     for (f in K1_list) out <- out + f(tvec, param)
#     out
#   }
#   K2fun <- function(tvec, param) {
#     d <- length(tvec)
#     out <- matrix(0, d, d) * param[1]
#     for (f in K2_list) out <- out + f(tvec, param)
#     out
#   }
#   K3opfun <- function(tvec, param, v1, v2, v3) {
#     total <- 0 * param[1]
#     for (f in K3_list) total <- total + f(tvec, param, v1, v2, v3)
#     total
#   }
#   K4opfun <- function(tvec, param, v1, v2, v3, v4) {
#     total <- 0 * param[1]
#     for (f in K4_list) total <- total + f(tvec, param, v1, v2, v3, v4)
#     total
#   }
#
#   tiltingfun <- function(tvec, param) {
#     total <- 0 * param[1]
#     for (f in tilting_list) total <- total + f(tvec, param)
#     total
#   }
#
#   ineqfun <- function(tvec, param) {
#     pieces <- lapply(ineq_list, function(f) f(tvec, param))
#     total_size <- sum(lengths(pieces))
#     out <- numeric(total_size) * param[1]
#     if (total_size == 0L) return(out)
#     idx <- 1L
#     for (p in pieces) {
#       lp <- length(p)
#       if (lp > 0L) {
#         out[idx:(idx + lp - 1L)] <- p
#         idx <- idx + lp
#       }
#     }
#     out
#   }
#
#   hist_pieces <- vapply(
#     cgf_list,
#     function(cg) paste(cg$call_history, collapse = " -> "),
#     character(1)
#   )
#   combined_history <- paste0("[", paste(hist_pieces, collapse = ", "), "]")
#
#   createCGF(
#     K = Kfun, K1 = K1fun, K2 = K2fun,
#     K3operator = K3opfun, K4operator = K4opfun,
#     tilting_exponent = tiltingfun,
#     ineq_constraint  = ineqfun,
#     op_name = c(combined_history, "sumOfIndependentCGF_slow"),
#     ...
#   )
# }
#
# sumOfIndependentCGF_slow <- function(cgf_list, iidReps = "any", block_size = NULL, ...) {
#   if (is.null(iidReps)) iidReps <- "any"
#   base <- .sumOfIndependentCGF_internal_slow(cgf_list, ...)
#   iidReplicatesCGF(cgf = base, iidReps = iidReps, block_size = block_size)
# }
#
#
# ## ------------------------------------------------------------
# ##    A "fast" version (should match your patched file)
# ##    If you've already updated the package file, you can skip
# ##    and just set sumOfIndependentCGF_fast <- sumOfIndependentCGF.
# ## ------------------------------------------------------------
# sumOfIndependentCGF_fast <- sumOfIndependentCGF
#
#
# ## ============================================================
# ## correctness + speed of factored operators
# ## ============================================================
#
# set.seed(1)
# d <- 30
#
# lambda1 <- function(theta) rep(exp(theta[1]), d)
# lambda2 <- function(theta) rep(exp(theta[2]), d)
#
# cg1 <- PoissonModelCGF(lambda = lambda1, iidReps = 1)
# cg2 <- PoissonModelCGF(lambda = lambda2, iidReps = 1)
#
# cgf_slow <- sumOfIndependentCGF_slow(list(cg1, cg2), iidReps = "any")
# cgf_fast <- sumOfIndependentCGF_fast(list(cg1, cg2), iidReps = "any")
#
# theta <- c(log(2), log(3))
# tvec  <- rnorm(d, sd = 0.1)
#
# ## Basic cumulant checks
# stopifnot(isTRUE(all.equal(cgf_slow$K(tvec, theta),  cgf_fast$K(tvec, theta),  tol = 1e-12)))
# stopifnot(isTRUE(all.equal(cgf_slow$K1(tvec, theta), cgf_fast$K1(tvec, theta), tol = 1e-12)))
# stopifnot(isTRUE(all.equal(cgf_slow$K2(tvec, theta), cgf_fast$K2(tvec, theta), tol = 1e-12)))
#
# ## Factored operator checks (call private methods by whitelist)
# Q <- crossprod(matrix(rnorm(d*d), d, d)) + diag(d)*0.1
# eig <- eigen(Q, symmetric = TRUE)
# A <- eig$vectors
# dvals <- eig$values
#
# K4slow <- cgf_slow$.get_private_method("K4operatorAABB_factored")
# K4fast <- cgf_fast$.get_private_method("K4operatorAABB_factored")
#
# ABBCslow <- cgf_slow$.get_private_method("K3K3operatorAABBCC_factored")
# ABBCfast <- cgf_fast$.get_private_method("K3K3operatorAABBCC_factored")
#
# ABCslow <- cgf_slow$.get_private_method("K3K3operatorABCABC_factored")
# ABCfast <- cgf_fast$.get_private_method("K3K3operatorABCABC_factored")
#
# r1 <- K4slow(tvec, theta, A, dvals, A, dvals)
# r2 <- K4fast(tvec, theta, A, dvals, A, dvals)
# cat("K4 AABB factored rel.err =", abs(r1-r2)/max(1,abs(r1)), "\n")
#
# r1 <- ABBCslow(tvec, theta, A, dvals, A, dvals, A, dvals)
# r2 <- ABBCfast(tvec, theta, A, dvals, A, dvals, A, dvals)
# cat("K3K3 AABBCC factored rel.err =", abs(r1-r2)/max(1,abs(r1)), "\n")
#
# system.time(r1 <- ABCslow(tvec, theta, A, dvals, A, dvals, A, dvals))
# system.time(r2 <- ABCfast(tvec, theta, A, dvals, A, dvals, A, dvals))
# cat("K3K3 ABCABC factored rel.err =", abs(r1-r2)/max(1,abs(r1)), "\n")
#
# ## Speed comparison: ABCABC is usually the expensive one (O(r^3))
# cat("\nTiming ABCABC_factored (slow vs fast):\n")
# system.time(replicate(3, ABCslow(tvec, theta, A, dvals, A, dvals, A, dvals)))
# system.time(replicate(3, ABCfast(tvec, theta, A, dvals, A, dvals, A, dvals)))
#
#
#
#
#
#
#
#
# ## ============================================================
# ## Multinomial (reduced dimension) + saddlepoint MLE
# ## ============================================================
# ## NOTE: Full multinomial has singular K2 (one linear constraint),
# ## so we drop the last category and work in dimension d-1.
#
# set.seed(2)
# d  <- 5
# B  <- 25
# N1 <- 40
# N2 <- 60
#
# theta_true <- c(0.3, -0.2, 0.1, -0.4)  # length d-1
#
# prob_odds <- function(theta) exp(c(theta, 0))  # odds; function normalizes internally
#
# cg_full1 <- MultinomialModelCGF(n = adaptor(fixed_param = N1), prob_vec = prob_odds, iidReps = 1)
# cg_full2 <- MultinomialModelCGF(n = adaptor(fixed_param = N2), prob_vec = prob_odds, iidReps = 1)
#
# A_sel <- cbind(diag(d-1), 0)  # (d-1) x d selection matrix
# cg_red1 <- linearlyMappedCGF(cg_full1, matrix_A = A_sel)
# cg_red2 <- linearlyMappedCGF(cg_full2, matrix_A = A_sel)
#
# ## Sum of independent multinomials (same p, sizes N1 and N2)
# cg_sum_one <- sumOfIndependentCGF_fast(list(cg_red1, cg_red2), iidReps = 1)
# cg_sum_B   <- iidReplicatesCGF(cg_sum_one, iidReps = B, block_size = d-1)
#
# ## Simulate data: equivalent to Multinomial(N1+N2, p)
# p_true <- prob_odds(theta_true); p_true <- p_true / sum(p_true)
# X_full <- rmultinom(B, size = N1 + N2, prob = p_true)  # d x B
# Y      <- X_full[1:(d-1), , drop = FALSE]              # reduced obs (d-1) x B
#
# ## Closed-form exact MLE for multinomial probabilities:
# p_hat <- rowSums(X_full) / sum(X_full)
# theta_hat_exact <- log(p_hat[1:(d-1)] / p_hat[d])
#
# ## Saddlepoint MLE:
# ## find.saddlepoint.MLE treats columns of observed.data as iid blocks and flattens column-wise.
# fit_mult <- find.saddlepoint.MLE(
#   observed.data  = Y,
#   cgf            = cg_sum_B,
#   starting.theta = rep(0, d-1),
#   method         = "two_step",
#   std.error = TRUE,
#   discrepancy = TRUE
# )
#
# cat("\nMultinomial example:\n")
# cat("theta_true      =", paste(round(theta_true, 4), collapse=" "), "\n")
# cat("theta_hat_exact =", paste(round(theta_hat_exact, 4), collapse=" "), "\n")
# cat("theta_hat_spa   =", paste(round(fit_mult$MLEs.theta, 4), collapse=" "), "\n")
# fit_mult$std.error
# fit_mult$discrepancy
#
# ## ============================================================
# ## Example C: Saddlepoint MLE for a nontrivial sum:
# ##           Y = Pois(lambda) + Binom(n,p)
# ##           Compare SPA-MLE vs exact convolution MLE
# ## ============================================================
#
# set.seed(3)
# B <- 200
# n_fixed <- 12
#
# lambda_true <- 3.5
# p_true      <- 0.35
# theta_true  <- c(log(lambda_true), qlogis(p_true))  # unconstrained parametrization
#
# ## Build CGFs (single observation)
# cg_pois <- PoissonModelCGF(lambda = function(th) exp(th[1]), iidReps = 1)
# cg_bin  <- BinomialModelCGF(n = adaptor(fixed_param = n_fixed),
#                             prob = function(th) plogis(th[2]), iidReps = 1)
#
# cg_sum <- sumOfIndependentCGF_fast(list(cg_pois, cg_bin), iidReps = B, block_size = 1)
#
# ## Simulate data
# y <- rpois(B, lambda_true) + rbinom(B, n_fixed, p_true)
#
# ## SPA MLE
# system.time(fit_sum <- find.saddlepoint.MLE(
#   observed.data  = y,
#   cgf            = cg_sum,
#   starting.theta = c(log(mean(y) + 1e-3), 0),
#   lb.theta = c(0.00001, -900),
#   discrepancy = TRUE,
#   method  = "two_step"
# ))
#
# system.time(fit_sum_slow <- find.saddlepoint.MLE(
#   observed.data  = y,
#   cgf            = sumOfIndependentCGF_slow(list(cg_pois, cg_bin), iidReps = B, block_size = 1),
#   starting.theta = c(log(mean(y) + 1e-3), 0),
#   lb.theta = c(0.00001, -900),
#   discrepancy = TRUE,
#   method = "constrained"
# ))
#
# ## Exact convolution log-likelihood for comparison
# negll_exact <- function(theta) {
#   lam <- exp(theta[1])
#   p   <- plogis(theta[2])
#
#   ll <- 0
#   for (yy in y) {
#     kmax <- min(n_fixed, yy)
#     ks <- 0:kmax
#     pm <- dbinom(ks, size = n_fixed, prob = p) * dpois(yy - ks, lambda = lam)
#     ll <- ll + log(sum(pm))
#   }
#   -ll
# }
#
# fit_exact <- nlminb(start = c(log(mean(y)), 0), objective = negll_exact)
#
# cat("\nPoisson + Binomial sum example:\n")
# cat("theta_true        =", paste(round(theta_true, 4), collapse=" "), "\n")
# cat("theta_hat_spa     =", paste(round(fit_sum$MLEs.theta, 4), collapse=" "), "\n")
# cat("theta_hat_spa_slow  =", paste(round(fit_sum_slow$MLEs.theta, 4), collapse=" "), "\n")
# cat("theta_hat_exactLL =", paste(round(fit_exact$par, 4), collapse=" "), "\n")
#
# fit_exact$par - fit_sum$MLEs.theta
# fit_sum$discrepancy
#
#
