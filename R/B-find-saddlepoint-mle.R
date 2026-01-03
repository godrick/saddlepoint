# R/find-saddlepoint-mle.R
# Objects: find.saddlepoint.MLE
#
#
# - Default behaviour is unchanged:
#     method = "constrained"
#   performs joint constrained optimisation over (tvec, theta) using nloptr.
#
# - New behavior:
#     method = "two_step"
#   optimizes over theta only and computes t-hat(theta) internally using the
#   RTMB/Newton pathway already used by compute.spa.negll / create_spa_taped_fun.
#
# Notes:
# - User inequality constraints (user.ineq.constraint.function) are assumed to be
#   constraints on theta only, matching get.ineq.constraint.function().
# - CGF domain constraints (cgf$ineq_constraint) are always enforced in
#   method="constrained". In method="two_step" they are enforced implicitly by
#   the internal t-hat solver; if the solver fails, the objective is penalised.


#' @noRd
configure.opts <- function(opts.user) {
  # Default optimizer options (nloptr / NLOPT_LD_SLSQP)
  opts.default <- list(
    algorithm = "NLOPT_LD_SLSQP",
    ftol_abs = 0,
    maxeval = 1e3,
    xtol_rel = 1.0e-12,
    print_level = 0
  )

  valid.option.names <- setdiff(names(opts.default), "algorithm") # exclude "algorithm"
  if (is.null(opts.user)) opts.user <- list()

  if (any(names(opts.user) == "")) {
    stop("All elements in opts.user must have names. Valid options are: ",
         paste(valid.option.names, collapse = ", "))
  }
  if (!all(names(opts.user) %in% valid.option.names)) {
    stop("Invalid option name(s) provided. Valid options are: ",
         paste(valid.option.names, collapse = ", "))
  }
  modifyList(opts.default, opts.user)
}


#' @noRd
.compute_t_hat_for_reporting <- function(observed.data,
                                        cgf,
                                        theta_hat,
                                        starting.tvec,
                                        lb.tvec,
                                        ub.tvec) {

  # Analytic saddlepoint (if provided by CGF)
  if (isTRUE(cgf$has_analytic_tvec_hat())) {
    t_hat <- as.numeric(cgf$analytic_tvec_hat(observed.data, theta_hat))
    return(list(t_hat = t_hat, t_source = "analytic"))
  }

  # Newton preferred
  t_hat <- tryCatch({
    G <- build_tvec_hat_newton_tape(
      cgf = cgf,
      y = observed.data,
      theta_init = theta_hat,
      t_init = starting.tvec
    )
    as.numeric(G(theta_hat))
  }, error = function(e) NULL)

  if (!is.null(t_hat) && all(is.finite(t_hat))) {
    # Guard against a Newton failure that returns a finite but inaccurate point
    res <- max(abs(cgf$K1(t_hat, theta_hat) - observed.data))
    if (is.finite(res) && res < 1e-4) return(list(t_hat = t_hat, t_source = "newton"))
  }


  # Fallback numeric solver
  sol <- saddlepoint.solve(
    theta = theta_hat,
    y = observed.data,
    cgf = cgf,
    starting.tvec = starting.tvec,
    lb = lb.tvec,
    ub = ub.tvec
  )
  # list(t_hat = as.numeric(sol$tvec_hat_vals), t_source = "solver")
  list(t_hat = as.numeric(sol), t_source = "solver")

}


#' @noRd
.fit_saddlepoint_MLE_constrained <- function(observed.data,
                                            cgf,
                                            starting.theta,
                                            lb.theta,
                                            ub.theta,
                                            starting.tvec,
                                            lb.tvec,
                                            ub.tvec,
                                            std.error,
                                            discrepancy,
                                            user.ineq.constraint.function,
                                            opts.user,
                                            zeroth.order) {

  # Combine tvec and theta arguments into a single vector for the optimizer
  initial_x <- c(starting.tvec, starting.theta)

  # Define objective, equality constraint and inequality constraint functions
  objective.function <- get.saddlepoint.nll.function(
    tvec = starting.tvec,
    theta = starting.theta,
    cgf = cgf
  )
  eq.constraint.function <- get.saddlepoint.eq.constraint.function(
    tvec = starting.tvec,
    theta = starting.theta,
    observed.data = observed.data,
    cgf = cgf
  )
  ineq.constraint.function <- get.ineq.constraint.function(
    tvec = starting.tvec,
    theta = starting.theta,
    cgf = cgf,
    user.ineq.constraint.function = user.ineq.constraint.function
  )

  # Override objective function if zeroth order method is selected
  if (isTRUE(zeroth.order)) {
    objective.function <- get.zeroth.saddlepoint.nll.function(
      tvec = starting.tvec,
      theta = starting.theta,
      cgf = cgf
    )
  }

  # Configure optimizer options
  opts <- configure.opts(opts.user)

  # Find the maximum likelihood estimates of the parameters using NLOPT
  MLEs <- nloptr::nloptr(
    x0          = initial_x,
    eval_f      = objective.function,
    eval_g_eq   = eq.constraint.function,
    eval_g_ineq = ineq.constraint.function,
    opts        = opts,
    lb          = c(lb.tvec, lb.theta),
    ub          = c(ub.tvec, ub.theta)
  )

  MLEs.tvec  <- head(MLEs$solution, length(lb.tvec))
  MLEs.theta <- tail(MLEs$solution, length(lb.theta))

  if (MLEs$status < 0) {
    cat("Estimates: ", MLEs.theta, "\n")
    stop("Optimization failed: ", MLEs$message)
  }

  # Calculate standard errors of MLEs
  if (isTRUE(std.error) || isTRUE(discrepancy)) {
    res <- compute.std.error(
      observed.data    = observed.data,
      estimated.tvec   = MLEs.tvec,
      estimated.theta  = MLEs.theta,
      cgf              = cgf,
      zeroth.order     = zeroth.order,
      non.saddlepoint.negll.function = NULL
    )
    MLEs$std.error <- res$std.error
    MLEs$inverse.hessian <- res$inverse.hessian
  }

  if (isTRUE(discrepancy)) {
    spa_method <- if (isTRUE(zeroth.order)) "zeroth" else "standard"
    out_ <- compute.saddlepointLL.correction(
      parameter_vector = MLEs.theta,
      observed.data    = observed.data,
      cgf              = cgf,
      tvec.hat         = MLEs.tvec,
      gradient         = TRUE,
      hessian          = FALSE,
      spa_method       = spa_method
    )
    discr <- MLEs$inverse.hessian %*% out_$gradient
    MLEs$discrepancy <- as.vector(discr)
  }

  # MLEs$MLEs.tvec  <- MLEs.tvec
  # MLEs$MLEs.theta <- MLEs.theta
  # MLEs$method <- "constrained"
  # MLEs

  MLEs$MLEs.tvec  <- MLEs.tvec
  MLEs$MLEs.theta <- MLEs.theta
  MLEs$method <- "constrained"
  MLEs$optimizer <- "nloptr"
  MLEs$outer_iterations <- if (!is.null(MLEs$iterations)) MLEs$iterations else NA_integer_
  MLEs$outer_iterations_type <- "nloptr major iterations (joint over tvec,theta)"
  # if (!is.null(MLEs$evaluations)) MLEs$outer_evaluations <- MLEs$evaluations
  MLEs
}



#' @noRd
.fit_saddlepoint_MLE_two_step <- function(observed.data,
                                         cgf,
                                         starting.theta,
                                         lb.theta,
                                         ub.theta,
                                         starting.tvec,
                                         lb.tvec,
                                         ub.tvec,
                                         std.error,
                                         discrepancy,
                                         user.ineq.constraint.function,
                                         opts.user,
                                         control.nlminb,
                                         zeroth.order) {

  has_user_ineq <- !is.null(user.ineq.constraint.function)

  # SPA method names for the RTMB tape builder:
  # create_spa_taped_fun() expects "negll_standard" / "negll_zeroth".
  spa_method_tape <- if (isTRUE(zeroth.order)) "negll_zeroth" else "negll_standard"


  # Detect whether the CGF imposes a domain constraint on tvec.
  # If so, the unconstrained RTMB/Newton path ("auto") may step outside the
  # domain; we default to the constrained numeric solver wrapped as an atomic.
  ineq_try <- tryCatch(cgf$ineq_constraint(starting.tvec, starting.theta),
                       error = function(e) e)
  if (inherits(ineq_try, "error")) {
    stop(
      "Unable to evaluate cgf$ineq_constraint(starting.tvec, starting.theta). ",
      "Maybe a dimension mismatch (block_size/iidReps) ??? .\n",
      "Original error: ", conditionMessage(ineq_try)
    )
  }

  has_cgf_ineq <- (length(ineq_try) > 0L)
  has_analytic_t_hat <- isTRUE(cgf$has_analytic_tvec_hat())
  # If an analytic t-hat exists, keep the fast auto/analytic pathway.
  # Otherwise, if the CGF is constrained, fall back to the constrained numeric solver.
  tvec_source_use <- if (has_analytic_t_hat) "auto" else if (has_cgf_ineq) "solver_atomic" else "auto"


  if (identical(tvec_source_use, "solver_atomic")) {
    if (isTRUE(has_user_ineq)) {
      warning(
        "method='two_step' was requested, but this CGF has inequality constraints on tvec. ",
        "In the current implementation, the two-step approach will be extremely slow. ",
        "For constrained CGFs, consider method='constrained' for speed/robustness. ",
        "Performance for this case may improve in future versions.",
        call. = FALSE, immediate. = TRUE
      )
    } else {
      warning(
        "method='two_step' was requested, but this CGF has inequality constraints on tvec. ",
        "In the current implementation, the two-step approach will be extremely slow. ",
        "For constrained CGFs, consider method='constrained' for speed/robustness. ",
        "Performance for this case may improve in future versions.",
        call. = FALSE, immediate. = TRUE
      )
    }
  }

  # Wrap the t-hat numeric solver so (optional) user-supplied tvec bounds are honoured.
  solver_use <- function(theta, y, cgf, starting.tvec) {
    saddlepoint.solve(theta = theta, y = y, cgf = cgf,
                      starting.tvec = starting.tvec,
                      lb = lb.tvec, ub = ub.tvec)
  }



  taped_spa <- create_spa_taped_fun(
    param_vec     = starting.theta,
    observed.data = observed.data,
    cgf           = cgf,
    spa_method    = spa_method_tape,
    tvec_source   = tvec_source_use,
    solver_fun    = solver_use,
    user_tvec     = NULL,
    newton_t_init = starting.tvec,
    gradient      = TRUE,
    hessian       = FALSE
  )



  cache_ <- new.env(parent = emptyenv())
  cache_$theta <- NULL
  cache_$res <- NULL

  .eval_cached <- function(theta) {
    theta <- as.numeric(theta)
    if (!is.null(cache_$theta) && isTRUE(all.equal(theta, cache_$theta, tolerance = 0))) {
      return(cache_$res)
    }
    r <- tryCatch(taped_spa(theta), error = function(e) NULL)
    if (is.null(r) || length(r$vals) != 1 || !is.finite(r$vals)) {
      r <- list(vals = 1e100, gradient = rep(1e-06, length(theta)))
    }
    cache_$theta <- theta
    cache_$res <- r
    r
  }




  # outer optimizer selection
  # If the user supplies inequality constraints, we use nloptr.
  # Otherwise, nlminb
  if (!has_user_ineq) {

    # If the user supplied a partial list, merge it into sensible defaults.
    if (is.null(control.nlminb)) control.nlminb <- list()
    ctrl_default <- list(eval.max = 200, iter.max = 150, rel.tol = 1.0e-10)
    nlminb_control <- modifyList(ctrl_default, control.nlminb)

    obj_fun  <- function(theta) as.numeric(.eval_cached(theta)$vals)
    grad_fun <- function(theta) as.numeric(.eval_cached(theta)$gradient)

    res_theta <- stats::nlminb(
      start     = starting.theta,
      objective = obj_fun,
      gradient  = grad_fun,
      lower     = lb.theta,
      upper     = ub.theta,
      control   = nlminb_control
    )

    theta_hat <- as.numeric(res_theta$par)
    opt_value <- as.numeric(res_theta$objective)


    # if (!isTRUE(res_theta$convergence == 0)) {
    #   cat("Estimates: ", theta_hat, "\n")
    #   stop("Optimization failed (nlminb): ", res_theta$message)
    # }

    # Compute t-hat at theta_hat for downstream functions
    t_hat_list <- .compute_t_hat_for_reporting(
      observed.data = observed.data,
      cgf = cgf,
      theta_hat = theta_hat,
      starting.tvec = starting.tvec,
      lb.tvec = lb.tvec,
      ub.tvec = ub.tvec
    )
    t_hat <- t_hat_list$t_hat

    # Mimic the nloptr output structure
    out <- list(
      solution  = c(t_hat, theta_hat),
      objective = opt_value,
      status    = res_theta$convergence,
      message   = res_theta$message
    )
    out$MLEs.tvec   <- t_hat
    out$MLEs.theta  <- theta_hat
    out$tvec.source <- t_hat_list$t_source
    out$method      <- "two_step"
    out$optimizer   <- "nlminb"



    out$outer_iterations <- res_theta$iterations
    out$outer_iterations_type <- "nlminb outer iterations"
    # if (!is.null(res_theta$evaluations)) out$outer_evaluations <- res_theta$evaluations
    out$cgf_constrained <- has_cgf_ineq
    out$has_user_ineq <- has_user_ineq
    out$tvec_source_tape <- tvec_source_use
    out$has_analytic_t_hat <- has_analytic_t_hat


    # Standard errors / discrepancy (same downstream code as constrained)
    if (isTRUE(std.error) || isTRUE(discrepancy)) {
      se <- compute.std.error(
        observed.data    = observed.data,
        estimated.tvec   = out$MLEs.tvec,
        estimated.theta  = out$MLEs.theta,
        cgf              = cgf,
        zeroth.order     = zeroth.order,
        non.saddlepoint.negll.function = NULL
      )
      out$std.error <- se$std.error
      out$inverse.hessian <- se$inverse.hessian
    }

    if (isTRUE(discrepancy)) {
      spa_method <- if (isTRUE(zeroth.order)) "zeroth" else "standard"
      corr <- compute.saddlepointLL.correction(
        parameter_vector = out$MLEs.theta,
        observed.data    = observed.data,
        cgf              = cgf,
        tvec.hat         = out$MLEs.tvec,
        gradient         = TRUE,
        hessian          = FALSE,
        spa_method       = spa_method
      )
      out$discrepancy <- as.vector(out$inverse.hessian %*% corr$gradient)
    }

    return(out)
  }

  # nloptr path (theta-only constraints
  # user.ineq.constraint.function(theta) must return list(constraints, jacobian),
  # with feasibility defined by constraints <= 0 (NLOPT convention).
  theta_obj <- function(theta) {
    r <- .eval_cached(theta)
    list(objective = as.numeric(r$vals), gradient = as.numeric(r$gradient))
  }

  opts <- configure.opts(opts.user)

  res_theta <- nloptr::nloptr(
    x0          = starting.theta,
    eval_f      = theta_obj,
    eval_g_ineq = user.ineq.constraint.function,
    lb          = lb.theta,
    ub          = ub.theta,
    opts        = opts
  )


  theta_hat <- as.numeric(res_theta$solution)

  if (res_theta$status < 0) {
    cat("Estimates: ", theta_hat, "\n")
    stop("Optimization failed (nloptr): ", res_theta$message)
  }

  # Compute t-hat at theta_hat
  t_hat_list <- .compute_t_hat_for_reporting(
    observed.data = observed.data,
    cgf = cgf,
    theta_hat = theta_hat,
    starting.tvec = starting.tvec,
    lb.tvec = lb.tvec,
    ub.tvec = ub.tvec
  )

  res_theta$MLEs.tvec   <- t_hat_list$t_hat
  res_theta$MLEs.theta  <- theta_hat
  res_theta$tvec.source <- t_hat_list$t_source
  res_theta$method      <- "two_step"
  res_theta$optimizer   <- "nloptr"

  res_theta$outer_iterations <- if (!is.null(res_theta$iterations)) res_theta$iterations else NA_integer_
  res_theta$outer_iterations_type <- "nloptr major iterations (theta-only)"
  if (!is.null(res_theta$evaluations)) res_theta$outer_evaluations <- res_theta$evaluations
  res_theta$cgf_constrained <- has_cgf_ineq
  res_theta$has_user_ineq <- has_user_ineq
  res_theta$tvec_source_tape <- tvec_source_use
  res_theta$has_analytic_t_hat <- has_analytic_t_hat




  if (isTRUE(std.error) || isTRUE(discrepancy)) {
    se <- compute.std.error(
      observed.data    = observed.data,
      estimated.tvec   = res_theta$MLEs.tvec,
      estimated.theta  = res_theta$MLEs.theta,
      cgf              = cgf,
      zeroth.order     = zeroth.order,
      non.saddlepoint.negll.function = NULL
    )
    res_theta$std.error <- se$std.error
    res_theta$inverse.hessian <- se$inverse.hessian
  }

  if (isTRUE(discrepancy)) {
    spa_method <- if (isTRUE(zeroth.order)) "zeroth" else "standard"
    corr <- compute.saddlepointLL.correction(
      parameter_vector = res_theta$MLEs.theta,
      observed.data    = observed.data,
      cgf              = cgf,
      tvec.hat         = res_theta$MLEs.tvec,
      gradient         = TRUE,
      hessian          = FALSE,
      spa_method       = spa_method
    )
    res_theta$discrepancy <- as.vector(res_theta$inverse.hessian %*% corr$gradient)
  }

  res_theta
}


#' @title Find maximum likelihood estimates using the saddlepoint likelihood.
#'
#' @description
#' This function finds maximum likelihood estimates (MLEs) under a saddlepoint
#' approximation, using constrained or two-step optimization.
#'
#' By default (\code{method="constrained"}) it performs the joint
#' constrained optimisation over \eqn{(t,\theta)}:
#' \deqn{ \hat t, \hat\theta = \arg\min_{(t,\theta)} \; \mathrm{NLL}_{\mathrm{SPA}}(t,\theta)
#' \quad \text{s.t.}\quad K_1(t;\theta) = y, }
#' using \code{NLOPT_LD_SLSQP}.
#'
#' Optionally, you may set \code{method="two_step"} to optimise over \eqn{\theta}
#' only, where \eqn{\hat t(\theta)} is computed internally for each \eqn{\theta}
#' using the RTMB/Newton pathway already implemented in \code{compute.spa.negll}.
#' This typically reduces the outer problem dimension and MAY be much faster.
#'
#' @details
#' - Observed Data: \code{observed.data} can be a numeric vector or a matrix/data frame.
#'   If it's a matrix/data frame, columns are treated as i.i.d. replicate blocks.
#'
#' - Zeroth-order: If \code{zeroth.order = TRUE}, the objective omits the
#'   saddlepoint correction term and uses the zeroth-order approximation.
#   \deqn{\hat t^T y - K_Y(\hat t;\theta).}
#'
#' - Inequality constraints:
#'   \itemize{
#'     \item Constraints from \code{cgf$ineq_constraint(tvec,theta)} are enforced in
#'       \code{method="constrained"}.
#'     \item In \code{method="two_step"}, \code{cgf$ineq_constraint} is enforced implicitly
#'       by the internal \eqn{\hat t(\theta)} solver; candidates where no feasible
#'       \eqn{\hat t(\theta)} can be found are penalized.
#'     \item \code{user.ineq.constraint.function}, if supplied, is interpreted as a
#'       constraint on \eqn{\theta} only. It must have signature
#'       \code{function(theta) -> list(constraints, jacobian)} with feasibility defined by
#'       \code{constraints <= 0} (NLOPT convention).
#'   }
#'
#' @param observed.data A numeric vector, matrix, or data frame containing the observed data.
#' @param cgf A CGF object corresponding to the distribution of the observed data.
#' @param starting.theta A numeric vector of starting values for model parameters.
#' @param lb.theta A numeric vector of lower bounds for \code{theta}. Defaults to \code{-Inf}.
#' @param ub.theta A numeric vector of upper bounds for \code{theta}. Defaults to \code{Inf}.
#' @param starting.tvec A numeric vector of starting values for the saddlepoint \eqn{t}.
#'   Defaults to \code{rep(0, length(observed.data))}.
#' @param lb.tvec Lower bounds for \eqn{t}. Defaults to \code{-Inf}.
#' @param ub.tvec Upper bounds for \eqn{t}. Defaults to \code{Inf}.
#' @param std.error Logical. If \code{TRUE}, compute standard errors of the MLEs. Defaults to \code{FALSE}.
#' @param discrepancy Logical. If \code{TRUE}, compute the discrepancy approximation. Defaults to \code{FALSE}.
#' @param user.ineq.constraint.function Optional user-defined inequality constraints on \eqn{\theta}.
#'   Must return a list with elements \code{constraints} and \code{jacobian}, with feasibility
#'   defined by \code{constraints <= 0}.
#' @param opts.user A named list of options for the \pkg{nloptr} optimizer.
#'   By default: \code{ftol_abs = 0}, \code{maxeval = 1e4}, \code{xtol_rel = 1e-7}, \code{print_level = 0}.
#'   The algorithm is fixed to \code{"NLOPT_LD_SLSQP"}.
#' @param zeroth.order Logical; if \code{TRUE} use the zeroth-order SPA objective. Defaults to \code{FALSE}.
#' @param method Optimisation strategy: \code{"constrained"} (default; backward compatible)
#'   or \code{"two_step"}.
#' @param control.nlminb Optional control list passed to \code{nlminb(..., control=...)} when
#'   \code{method="two_step"} and there are no user inequality constraints. If \code{NULL},
#'   defaults are \code{eval.max=200}, \code{iter.max=150}, \code{rel.tol=1e-10}.
#'
#' @return A list containing (at least):
#'   \itemize{
#'     \item \code{MLEs.tvec, MLEs.theta}: the estimated saddlepoint \eqn{\hat{t}} and parameter vector \eqn{\hat{\theta}}.
#'     \item \code{std.error}, \code{inverse.hessian}: if \code{std.error=TRUE}.
#'     \item \code{discrepancy}: if \code{discrepancy=TRUE}.
#'     \item \code{solution, status, message}: raw optimiser output (structure depends on optimiser).
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1); x <- rgamma(50, shape = 10, rate = 0.5)
#' res <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF,
#'                             starting.theta = c(1,1), std.error=TRUE)
#' res$MLEs.theta
#' res$std.error
#'
#' ## Two-step optimization:
#' res2 <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF,
#'                              starting.theta = c(1,1), method="two_step", std.error=TRUE)
#' res2$MLEs.theta
#' res2$std.error
#' }
#' @export
find.saddlepoint.MLE <- function(observed.data,
                                 cgf,
                                 starting.theta,
                                 lb.theta = rep(-Inf, times = length(starting.theta)),
                                 ub.theta = rep( Inf, times = length(starting.theta)),
                                 starting.tvec = rep(0, times = length(as.numeric(observed.data))),
                                 lb.tvec = rep(-Inf, times = length(as.numeric(observed.data))),
                                 ub.tvec = rep( Inf, times = length(as.numeric(observed.data))),
                                 std.error = FALSE,
                                 discrepancy = FALSE,
                                 user.ineq.constraint.function = NULL,
                                 opts.user = list(ftol_abs = 0, maxeval = 1e4, xtol_rel = 1.0e-7, print_level = 0),
                                 zeroth.order = FALSE,
                                 method = c("constrained", "two_step"),
                                 control.nlminb = NULL) {

  method <- match.arg(method)

  # ------------------------------- data handling -------------------------------
  if (is.matrix(observed.data) || is.data.frame(observed.data)) {
    message("Treating columns of 'observed.data' as i.i.d. replicate blocks.")
    observed.data <- as.numeric(observed.data)
  }
  observed.data <- as.numeric(observed.data)

  # ------------------------------- checks ------------------------------------
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  if (!is.numeric(starting.theta)) stop("starting.theta must be numeric")
  if (!is.numeric(lb.theta) || !is.numeric(ub.theta)) stop("lb.theta and ub.theta must be numeric")
  if (!is.numeric(starting.tvec) || !is.numeric(lb.tvec) || !is.numeric(ub.tvec)) {
    stop("starting.tvec, lb.tvec, ub.tvec must be numeric")
  }

  if (length(lb.theta) != length(starting.theta) || length(ub.theta) != length(starting.theta)) {
    stop("lb.theta or ub.theta has an incorrect length")
  }
  if (any(starting.theta < lb.theta) || any(starting.theta > ub.theta)) {
    stop("starting.theta is not within the bounds specified")
  }

  if (length(starting.tvec) != length(observed.data)) {
    stop("Size of observed.data and starting.tvec arguments do not match")
  }
  if (length(lb.tvec) != length(starting.tvec) || length(ub.tvec) != length(starting.tvec)) {
    stop("lb.tvec or ub.tvec has an incorrect length")
  }

  # ------------------------------- dispatch -----------------------------------
  if (identical(method, "constrained")) {
    return(.fit_saddlepoint_MLE_constrained(
      observed.data = observed.data,
      cgf = cgf,
      starting.theta = starting.theta,
      lb.theta = lb.theta,
      ub.theta = ub.theta,
      starting.tvec = starting.tvec,
      lb.tvec = lb.tvec,
      ub.tvec = ub.tvec,
      std.error = std.error,
      discrepancy = discrepancy,
      user.ineq.constraint.function = user.ineq.constraint.function,
      opts.user = opts.user,
      zeroth.order = zeroth.order
    ))
  }

  .fit_saddlepoint_MLE_two_step(
    observed.data = observed.data,
    cgf = cgf,
    starting.theta = starting.theta,
    lb.theta = lb.theta,
    ub.theta = ub.theta,
    starting.tvec = starting.tvec,
    lb.tvec = lb.tvec,
    ub.tvec = ub.tvec,
    std.error = std.error,
    discrepancy = discrepancy,
    user.ineq.constraint.function = user.ineq.constraint.function,
    opts.user = opts.user,
    control.nlminb = control.nlminb,
    zeroth.order = zeroth.order
  )
}
