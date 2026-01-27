# R/compute-spa-negll.R
# t-hat and SPA likelihood with multiple tvec sources (user/analytic/newton/solver_atomic)


# get_nonAD_tvec_hat_vals <- function(parameter_vector,
#                                     observed.data,
#                                     cgf,
#                                     user_tvec = NULL,
#                                     tvec_source = c("auto","user","analytic","newton","solver_atomic"),
#                                     solver_fun = saddlepoint.solve,
#                                     ...) {
#   tvec_source <- match.arg(tvec_source)
#   if (!is.null(user_tvec) || tvec_source == "user") return(user_tvec)
#   if (isTRUE(cgf$has_analytic_tvec_hat) || tvec_source == "analytic")
#     return(cgf$analytic_tvec_hat(observed.data, parameter_vector))
#   # # Otherwise plain numeric solve (no AD)
#   # solver_fun(theta = parameter_vector, y = observed.data, cgf = cgf, ...)
#
#
#   tvec_hat_vals <- tryCatch({
#     # robust: soft barrier by default for a single evaluation
#     Gtmp <- build_tvec_hat_newton_tape(
#       cgf = cgf, y = observed.data, theta_init = parameter_vector,
#       t_init = if (is.null(tvec.hat)) rep(0, length(observed.data)) else tvec.hat,
#       barrier = "soft"
#     )
#     Gtmp(parameter_vector)
#   }, error = function(e) {
#     # fallback to SLSQP if RTMB Newton fails for some reason
#     message("Falling back to saddlepoint.solve(): ", conditionMessage(e))
#     saddlepoint.solve(theta = parameter_vector, y = observed.data, cgf = cgf, ...)
#   })
#   tvec_hat_vals
#
# }




# -----------------------------------------------------------------------------
# Fast, no-AD path to get t-hat values
# -----------------------------------------------------------------------------
#' @noRd
get_nonAD_tvec_hat_vals <- function(parameter_vector,
	                                    observed.data,
	                                    cgf,
	                                    user_tvec = NULL,
	                                    tvec_source = c("auto","user","analytic","newton","solver_atomic"),
	                                    solver_fun = saddlepoint.solve,
	                                    ...) {
	  tvec_source <- match.arg(tvec_source)
	  if (tvec_source == "user") {
	    if (is.null(user_tvec)) stop("tvec_source='user' requires 'user_tvec'.", call. = FALSE)
	    return(user_tvec)
	  }
	  if (!is.null(user_tvec)) return(user_tvec)

	  # analytic t-hat (if available, or explicitly requested)
	  if (tvec_source == "analytic") {
	    if (!isTRUE(cgf$has_analytic_tvec_hat) || is.null(cgf$analytic_tvec_hat)) {
	      stop("tvec_source='analytic' requires this CGF to provide 'analytic_tvec_hat'.", call. = FALSE)
	    }
	    return(cgf$analytic_tvec_hat(observed.data, parameter_vector))
	  }
	  if (isTRUE(cgf$has_analytic_tvec_hat) && !is.null(cgf$analytic_tvec_hat)) {
	    return(cgf$analytic_tvec_hat(observed.data, parameter_vector))
	  }


  # Pure numeric solve (no tapes). We prefer not to build a Newton tape here.
  # When the caller wants gradients/Hessians they should take the AD path.
  # For robustness this is the default for auto/newton/solver_atomic in *non-AD* context.
  solver_fun(theta = parameter_vector, y = observed.data, cgf = cgf, ...)
}



# -----------------------------------------------------------------------------
# Which SPA function do we evaluate?
# -----------------------------------------------------------------------------
#' @noRd
choose_spa_function <- function(spa_method, cgf) {
  allowed_methods <- c(
    "negll_standard",        # full first-order neg log-lik
    "negll_zeroth",          # zeroth-order neg log-lik
    "correction_standard",   # first-order correction term T
    "correction_zeroth"      # zeroth-order correction term -0.5 log det K2
  )
  if (!spa_method %in% allowed_methods) {
    stop(sprintf(
      "Unknown saddlepoint method '%s'. Allowed methods: %s",
      spa_method, paste(allowed_methods, collapse = ", ")
    ))
  }

  neg_ll_fun  <- cgf$.private_api$neg_ll          # function(tvec, theta)
  tilt_fun    <- cgf$.private_api$tilting_exponent
  funcT_first <- cgf$.private_api$func_T

  if (spa_method == "negll_standard") {
    return(neg_ll_fun)
  }
  if (spa_method == "negll_zeroth") {
    return( function(tvec, theta) -tilt_fun(tvec, theta) )
  }
  if (spa_method == "correction_standard") {
    return(funcT_first)
  }

  function(tvec, theta) {
      # K2_val <- cgf$K2(tvec, theta)
      # -0.5 * determinant(K2_val, logarithm = TRUE)$modulus
      # # -0.5 * my_logdet(as.matrix(K2_val))
    -0.5 * cgf$logdetK2(tvec, theta)
  }
}






# -----------------------------------------------------------------------------
# Build an AD-taped function for SPA using a chosen tvec source
# -----------------------------------------------------------------------------
#' @noRd
create_spa_taped_fun <- function(param_vec,
                                 observed.data,
                                 cgf,
                                 spa_method,     # "negll_standard", "negll_zeroth", "correction_*"
                                 tvec_source = c("auto","user","analytic","newton","solver_atomic"),
                                 user_tvec = NULL,
                                 gradient = FALSE,
                                 hessian = FALSE,
                                 solver_fun = saddlepoint.solve, # numeric solver if 'solver_atomic'
                                 newton_t_init = NULL,           # start for t in newton
                                 ...
) {
  tvec_source <- match.arg(tvec_source)
  if (hessian && !gradient) gradient <- TRUE
  chosen_spa_fn <- choose_spa_function(spa_method = spa_method, cgf = cgf)

  # User supplied tvec (one-time tape)
  if (tvec_source == "user" || (!is.null(user_tvec) && tvec_source == "auto")) {
    if (is.null(user_tvec)) stop("tvec_source='user' requires 'user_tvec'.")
    t_atomic <- make_user_tvec_atomic(user_tvec, cgf, theta_init = param_vec)

    local_fn <- function(par) {
      tvec_hat_vals <- t_atomic(par)
      chosen_spa_fn(tvec_hat_vals, par)
    }
    tape_obj <- RTMB::MakeTape(local_fn, x = param_vec)
    if (hessian) jacfun_obj <- tape_obj$jacfun()
    return(function(theta) {
      list(
        vals     = tape_obj(theta),
        gradient = if (gradient) as.vector(tape_obj$jacobian(theta)) else NULL,
        hessian  = if (hessian) jacfun_obj$jacobian(theta) else NULL
      )
    })
  }


  # Analytic t-hat
  if (tvec_source == "analytic") {
    if (!isTRUE(cgf$has_analytic_tvec_hat) || is.null(cgf$analytic_tvec_hat)) {
      stop("tvec_source='analytic' requires this CGF to provide 'analytic_tvec_hat'.", call. = FALSE)
    }
  }
  if (tvec_source == "analytic" || (isTRUE(cgf$has_analytic_tvec_hat) && tvec_source == "auto")) {
    local_fn <- function(par) {
      tvec_vals <- cgf$analytic_tvec_hat(observed.data, par)
      chosen_spa_fn(tvec_vals, par)
    }
    tape_obj <- RTMB::MakeTape(local_fn, x = param_vec)
    if (hessian) jacfun_obj <- tape_obj$jacfun()
    return(function(theta) {
      list(
        vals     = tape_obj(theta),
        gradient = if (gradient) as.vector(tape_obj$jacobian(theta)) else NULL,
        hessian  = if (hessian) jacfun_obj$jacobian(theta) else NULL
      )
    })
  }


  # Newton (recommended; fully AD-safe, higher-order)
  if (tvec_source == "newton" || tvec_source == "auto") {
    G_theta_to_t <- build_tvec_hat_newton_tape(
      cgf         = cgf,
      y           = observed.data,
      theta_init  = param_vec,
      t_init      = newton_t_init
    )
    local_fn <- function(par) {
      tvec_hat_vals <- G_theta_to_t(par)  # evaluated inside tape for grad_theta
      chosen_spa_fn(tvec_hat_vals, par)
    }
    tape_obj <- RTMB::MakeTape(local_fn, x = param_vec)
    if (hessian) jacfun_obj <- tape_obj$jacfun()
    return(function(theta) {
      list(
        vals     = tape_obj(theta),
        gradient = if (gradient) as.vector(tape_obj$jacobian(theta)) else NULL,
        hessian  = if (hessian) jacfun_obj$jacobian(theta) else NULL
      )
    })
  }




  # External numeric solver wrapped as an atomic (AD-safe)
  # This gives grad_theta through the implicit-function formula
  if (tvec_source == "solver_atomic") {
    t_atomic <- make_solver_tvec_atomic(
      cgf        = cgf,
      y          = observed.data,
      solver_fun = solver_fun,
      theta_init = param_vec,
      t_init     = newton_t_init
    )
    local_fn <- function(par) {
      tvec_hat_vals <- t_atomic(par)
      chosen_spa_fn(tvec_hat_vals, par)
    }
    tape_obj <- RTMB::MakeTape(local_fn, x = param_vec)
    if (hessian) jacfun_obj <- tape_obj$jacfun()
    return(function(theta) {
      list(
        vals     = local_fn(theta),
        gradient = if (gradient) as.vector(tape_obj$jacobian(theta)) else NULL,
        hessian  = if (hessian) jacfun_obj$jacobian(theta) else NULL
      )
    })
  }

  stop("Unhandled tvec_source.")
}

# -----------------------------------------------------------------------------
# Public entry: compute the SPA neg log-lik (or zeroth-order), with derivatives
# -----------------------------------------------------------------------------
#' Compute the saddlepoint negative log-likelihood
#'
#' Computes either the "standard" or "zeroth-order" saddlepoint approximation to the
#' negative log–likelihood, with options for how the saddlepoint t-vector is obtained.
#'
#' @param parameter_vector Numeric vector of parameters.
#' @param observed.data Numeric vector of observations.
#' @param cgf A `CGF` object.
#' @param tvec.hat Optional numeric t-vector (one-time tape).
#' @param gradient,hessian Logical flags for derivatives.
#' @param spa_method "standard" or "zeroth".
#' @param tvec_source One of "auto","user","analytic","newton","solver_atomic".
#' @param solver_fun Numeric solver used when tvec_source="solver_atomic". Default `saddlepoint.solve`.
#' @param newton_t_init Optional start for t in the Newton solver path.
#' @param ... Passed to solver_fun (when used).
#' @return list(vals=..., gradient=..., hessian=...)
#' @export
compute.spa.negll <- function(parameter_vector,
                              observed.data,
                              cgf,
                              tvec.hat   = NULL,
                              gradient   = FALSE,
                              hessian    = FALSE,
                              spa_method = "standard",
                              tvec_source = c("auto","user","analytic","newton","solver_atomic"),
                              solver_fun  = saddlepoint.solve,
                              newton_t_init = NULL,
                              ...) {
  if (!inherits(cgf, "CGF")) stop("`cgf` must be an object of class CGF")
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data))    stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be logical(1).")
  if (!is.logical(hessian)  || length(hessian)  != 1) stop("`hessian` must be logical(1).")

  # Normalize method labels
  if (spa_method == "standard") spa_method <- "negll_standard"
  else if (spa_method == "zeroth") spa_method <- "negll_zeroth"

  tvec_source <- match.arg(tvec_source)

  # If no derivatives requested -> direct numeric path fastest
  if (!gradient && !hessian) {
    tvec_hat_vals <- get_nonAD_tvec_hat_vals(
      parameter_vector = parameter_vector,
      observed.data    = observed.data,
      cgf              = cgf,
      user_tvec        = tvec.hat,
      tvec_source      = tvec_source,
      solver_fun       = solver_fun,
      ...
    )
    chosen <- choose_spa_function(spa_method, cgf)
    val <- chosen(tvec_hat_vals, parameter_vector)[1]
    return(list(vals = val, gradient = NULL, hessian = NULL))
  }

  # Build tape suitable for chosen tvec source
  taped_fun <- create_spa_taped_fun(
    param_vec     = parameter_vector,
    observed.data = observed.data,
    cgf           = cgf,
    spa_method    = spa_method,
    tvec_source   = if (!is.null(tvec.hat) && tvec_source=="auto") "user" else tvec_source,
    user_tvec     = tvec.hat,
    gradient      = gradient,
    hessian       = hessian,
    solver_fun    = solver_fun,
    newton_t_init = newton_t_init,
    ...
  )
  taped_fun(parameter_vector)

}
