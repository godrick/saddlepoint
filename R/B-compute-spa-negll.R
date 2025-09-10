# R/compute-spa-negll.R
# t-hat and SPA likelihood with multiple tvec sources (user/analytic/newton/solver_atomic)

# -----------------------------------------------------------------------------
# Compute tvec values (no AD tape needed - plain numeric)
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
  if (!is.null(user_tvec) || tvec_source == "user") return(user_tvec)
  if (cgf$has_analytic_tvec_hat() || tvec_source == "analytic")
    return(cgf$analytic_tvec_hat(observed.data, parameter_vector))
  # Otherwise plain numeric solve (no AD)
  solver_fun(theta = parameter_vector, y = observed.data, cgf = cgf, ...)
}

# -----------------------------------------------------------------------------
# Which SPA function do we evaluate? (as in your original)
# -----------------------------------------------------------------------------
#' @noRd
choose_spa_function <- function(spa_method, cgf) {
  allowed_methods <- c(
    "negll_standard",  
    "negll_zeroth",    
    "correction_standard",  
    "correction_zeroth"
  )
  if (!spa_method %in% allowed_methods) {
    stop(sprintf(
      "Unknown saddlepoint method '%s'. Allowed methods: %s",
      spa_method, paste(allowed_methods, collapse = ", ")
    ))
  }
  
  neg_ll_fun  <- cgf$.get_private_method("neg_ll")          # function(tvec, theta)
  tilt_fun    <- cgf$.get_private_method("tilting_exponent")
  funcT_first <- cgf$.get_private_method("func_T")
  
  if (spa_method == "negll_standard") {
    neg_ll_fun
  } else if (spa_method == "negll_zeroth") {
    function(tvec, theta) -tilt_fun(tvec, theta)
  } else if (spa_method == "correction_standard") {
    funcT_first
  } else {
    function(tvec, theta) {
      K2_val <- cgf$K2(tvec, theta)
      -0.5 * determinant(K2_val, logarithm = TRUE)$modulus
    }
  }
}

# -----------------------------------------------------------------------------
# Build an AD-taped function for SPA using a chosen tvec source
# -----------------------------------------------------------------------------
#' @noRd
create_spa_taped_fun <- function(param_vec,
                                 observed.data,
                                 cgf,
                                 spa_method,     # "negll_standard", "negll_zeroth", etc
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
  
  # ---- User supplied tvec (one-time tape) -----------------------------------
  if (tvec_source == "user" || (!is.null(user_tvec) && tvec_source == "auto")) {
    if (is.null(user_tvec)) stop("tvec_source='user' requires 'user_tvec'.")
    # NOTE: pass 'theta_init = param_vec'
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
  
  
  
  # ---- Analytic t-hat --------------------------------------------------------
  if (tvec_source == "analytic" || (cgf$has_analytic_tvec_hat() && tvec_source == "auto")) {
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
  
  # ---- Newton (recommended; fully AD-safe, higher-order) ---------------------
  if (tvec_source == "newton" || tvec_source == "auto") {
    # Build once: a Tape G s.t. G(theta) = t_hat(theta), using Tape$newton. (RTMB manual) 
    G_theta_to_t <- build_tvec_hat_newton_tape(
      cgf         = cgf,
      y           = observed.data,
      theta_init  = param_vec,
      t_init      = newton_t_init
    )
    
    local_fn <- function(par) {
      tvec_hat_vals <- G_theta_to_t(par)  # *inside the tape* so derivatives propagate
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
  
  # ---- External numeric solver wrapped in a (Hessian-safe) atomic --------------
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
        vals     = tape_obj(theta),
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
