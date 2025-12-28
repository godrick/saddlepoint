


.rtmb_value_real <- function(x) {
  if (inherits(x, "advector")) x <- RTMB:::getValues(x)
  if (is.complex(x)) x <- Re(x)
  # if (is.complex(x)) stop("Unexpected complex value. Maybe AD context leakage ???....")
  x
  # as.numeric(x)
}



# # Strictly-feasible, damped Newton solver for the saddlepoint equation:
# #   K1(tvec; theta) = y
# # subject to:
# #   cgf$ineq_constraint(tvec, theta) <= 0
# # plus optional box bounds lb/ub.
# #
# # Globalization uses merit function psi(t)=0.5*||K1(t)-y||^2 (root-finding merit),
# #' @noRd
# .saddlepoint.newton.solve <- function(theta, y, cgf,
#                                       starting.tvec = rep(0, length(y)),
#                                       lb = rep(-Inf, length(y)),
#                                       ub = rep( Inf, length(y)),
#                                       tol = 1e-10,
#                                       maxit = 80,
#                                       interior_margin = 1e-10,
#                                       armijo_c1 = 1e-4,
#                                       step_fraction = 0.99,
#                                       step_shrink = 0.5,
#                                       max_backtrack = 60,
#                                       reg_init = 0,
#                                       reg_grow = 10,
#                                       reg_max = 1e12,
#                                       warn_residual = TRUE,
#                                       warn_tol = 1e-4,
#                                       verbose = FALSE,
#                                       return_info = FALSE) {
#
#   stopifnot(inherits(cgf, "CGF"))
#
#   theta <- as.numeric(.rtmb_value_real(theta))
#   y     <- as.numeric(.rtmb_value_real(y))
#   m     <- length(y)
#
#   tvec <- as.numeric(.rtmb_value_real(starting.tvec))
#   if (length(tvec) != m) stop("starting.tvec must have same length as y.")
#
#   lb <- as.numeric(lb); ub <- as.numeric(ub)
#   if (length(lb) != m || length(ub) != m) stop("lb/ub must have length(y).")
#
#   # Clamp to bounds (non-strict)
#   tvec <- pmin(pmax(tvec, lb), ub)
#
#   # Inequality constraints from CGF
#   g_fun <- function(t) {
#     .rtmb_value_real(cgf$ineq_constraint(t, theta))
#   }
#
#   feasible <- function(t) {
#     # Strictly interior wrt bounds if finite
#     if (any(is.finite(lb) & (t <= lb + interior_margin))) return(FALSE)
#     if (any(is.finite(ub) & (t >= ub - interior_margin))) return(FALSE)
#
#     g <- g_fun(t)
#     if (!length(g)) return(TRUE)
#     if (!all(is.finite(g))) return(FALSE)
#     max(g) < -interior_margin
#   }
#
#   # Ensure a strictly feasible start by shrinking toward a safe base point.
#   if (!feasible(tvec)) {
#     base <- pmin(pmax(rep(0, m), lb + interior_margin), ub - interior_margin)
#     alpha <- 1.0
#     ok <- FALSE
#     for (k in seq_len(80)) {
#       t_try <- (1 - alpha) * base + alpha * tvec
#       if (feasible(t_try)) { tvec <- t_try; ok <- TRUE; break }
#       alpha <- alpha * 0.5
#     }
#     if (!ok) stop("Could not find a strictly feasible starting tvec.")
#   }
#
#   converged <- FALSE
#   last_res  <- NA_real_
#   it_used   <- 0L
#
#   for (it in seq_len(maxit)) {
#     it_used <- it
#
#     F_ <- as.numeric(.rtmb_value_real(cgf$K1(tvec, theta))) - y
#     res <- max(abs(F_))
#     last_res <- res
#
#     if (!is.finite(res)) stop("Non-finite residual encountered.")
#     if (verbose) {
#       g <- g_fun(tvec)
#       gmax <- if (length(g)) max(g) else -Inf
#       cat(sprintf("it=%d  res=%.3e  max_g=%.3e\n", it, res, gmax))
#     }
#
#     if (res <= tol) { converged <- TRUE; break }
#
#     psi0    <- 0.5 * sum(F_ * F_)
#     Fnorm2  <- sum(F_ * F_)
#
#     step_accepted <- FALSE
#     reg <- reg_init
#
#     # Try increasing regularization if needed
#     for (reg_try in seq_len(8)) {
#
#       # Compute direction d
#       d <- tryCatch({
#         if (reg == 0) {
#           # Use K2_solve directly
#           -as.numeric(.rtmb_value_real(cgf$K2_solve(tvec, theta, F_)))
#         } else {
#           # Regularized fallback: build K2 matrix numerically and solve
#           H <- as.matrix(.rtmb_value_real(cgf$K2(tvec, theta)))
#           H <- H + reg * base::diag(nrow(H))
#           -as.numeric(base::solve(H, F_))
#         }
#       }, error = function(e) NULL)
#
#       if (is.null(d) || length(d) != m || any(!is.finite(d))) {
#         reg <- if (reg == 0) 1e-8 else reg * reg_grow
#         if (reg > reg_max) break
#         next
#       }
#
#       # Compute an initial alpha limited by bounds (keep strict interior)
#       alpha <- 1.0
#
#       idxU <- which(is.finite(ub) & d > 0)
#       if (length(idxU)) {
#         alpha <- min(alpha, min((ub[idxU] - interior_margin - tvec[idxU]) / d[idxU]))
#       }
#
#       idxL <- which(is.finite(lb) & d < 0)
#       if (length(idxL)) {
#         alpha <- min(alpha, min((lb[idxL] + interior_margin - tvec[idxL]) / d[idxL]))
#       }
#
#       alpha <- max(0, min(1, alpha))
#       alpha <- step_fraction * alpha
#
#       if (!(alpha > 0)) {
#         reg <- if (reg == 0) 1e-8 else reg * reg_grow
#         if (reg > reg_max) break
#         next
#       }
#
#       # Line search: enforce feasibility AND Armijo decrease on psi(t)=0.5||F_||^2
#       for (bt in seq_len(max_backtrack)) {
#         t_try <- tvec + alpha * d
#
#         if (!feasible(t_try)) {
#           alpha <- alpha * step_shrink
#           next
#         }
#
#         F_try <- as.numeric(.rtmb_value_real(cgf$K1(t_try, theta))) - y
#         psi_try <- 0.5 * sum(F_try * F_try)
#
#         # Armijo: psi(t+αd) <= psi(t) - c1 * α * ||F_||^2
#         if (is.finite(psi_try) && psi_try <= psi0 - armijo_c1 * alpha * Fnorm2) {
#           tvec <- t_try
#           step_accepted <- TRUE
#           break
#         }
#
#         alpha <- alpha * step_shrink
#         if (!(alpha > 0)) break
#       }
#
#       if (step_accepted) break
#
#       # If line search couldn't accept, increase regularization and retry direction
#       reg <- if (reg == 0) 1e-8 else reg * reg_grow
#       if (reg > reg_max) break
#     }
#
#     if (!step_accepted) {
#       stop(sprintf("Feasible Newton failed at it=%d (res=%.3e).", it, res))
#     }
#   }
#
#   if (!converged) {
#     stop(sprintf("Feasible Newton did not converge in %d iterations (final res=%.3e).",
#                  maxit, last_res))
#   }
#
#   if (warn_residual) {
#     rfinal <- max(abs(as.numeric(.rtmb_value_real(cgf$K1(tvec, theta))) - y))
#     if (is.finite(rfinal) && rfinal > warn_tol) {
#       warning(sprintf("Saddlepoint residual %.3e exceeds warn_tol=%.3e.", rfinal, warn_tol))
#     }
#   }
#
#   if (return_info) {
#     return(list(t_hat = tvec, converged = TRUE, iters = it_used, residual = last_res))
#   }
#   tvec
# }




