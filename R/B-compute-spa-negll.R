# R/compute-spa-negll.R
# Objects: 




#-----------------------------------------------------------------
# Maybe AD wrapper around `tvec_hat_for_ad()`
#
# This function ensures that if the `theta` being passed is an
# advector (for AD), it re-runs the solution for tvec with `tvec_hat_for_ad`,
# capturing the AD dependency. Otherwise, it just returns the numeric tvec.
#-----------------------------------------------------------------

#' Compute or Return a t-vector Depending on AD
#'
#' @description
#' Returns the values of the saddlepoint equation after capturing the AD dependency.
#'
#' @param theta A numeric vector or an AD vector representing parameters.
#' @param tvec Numeric vector. 
#' @param observed.data Numeric vector of observed data.
#' @param K1_fn A function giving the first derivative of the CGF.
#' @param K2_solve_fn A function. See `create_tvec_hat_K2_solve_fn()`.
#'
#' @return A vector (numeric or AD) representing the relevant `tvec`.
#' 
#' @noRd
tvec_hat <- function(
    theta,
    current_tvec,
    observed.data,
    K1_fn,
    K2_solve_fn
) {
  if (inherits(theta, "advector")) {
    # re-compute for AD
    tvec_hat_for_ad(
      theta        = theta,
      tvec         = current_tvec, 
      observations = observed.data,
      K1_fn        = K1_fn,
      K2_solve_fn = K2_solve_fn
    )
  } else {
    current_tvec
  }
}


#' @noRd
get_tvec_hat_vals <- function(parameter_vector,
                              observed.data,
                              cgf,
                              user_tvec = NULL,
                              ...) {
  if (!is.null(user_tvec)) return(user_tvec)
  if (cgf$has_analytic_tvec_hat()) return(cgf$analytic_tvec_hat(observed.data, parameter_vector))
  # Otherwise, solve numerically
  saddlepoint.solve(theta = parameter_vector,
                    y     = observed.data,
                    cgf   = cgf,
                    ...)
}


#' @noRd
choose_negll_fun <- function(method, cgf) {
  # “allowed” methods up front
  allowed_methods <- c("standard", "zeroth")  # can expand in future?????
  
  if (!method %in% allowed_methods) {
    stop(sprintf("Unknown saddlepoint method '%s'. Allowed methods: %s",
                 method, paste(allowed_methods, collapse=", ")))
  }
  
  neg_ll_fun   <- cgf$.get_private_method("neg_ll")
  tilting_fun  <- cgf$.get_private_method("tilting_exponent")
  
  if (method == "standard") {
    return(neg_ll_fun)
  } else {
    # method == "zeroth"
    # Return - tilting exponent
    return( function(tvec, parameter_vector) {-tilting_fun(tvec, parameter_vector)} )
  }
}




#' @noRd
get_spa_taped_fun <- function(param_vec,
                              observed.data,
                              cgf,
                              method = "standard",
                              user_tvec = NULL,
                              gradient = FALSE,
                              hessian = FALSE,
                              ... # additional arguments to saddlepoint.solve
) {
  if (hessian && !gradient) gradient <- TRUE
  chosen_negll <- choose_negll_fun(method = method, cgf = cgf)
  K2_solve_fn <- create_tvec_hat_K2_solve_fn(cgf)
  
  # The function that will be taped. Each time it's called with a parameter vector `par`,
  # it figures out how to compute tvec (user-supplied, analytic, or numeric).
  local_fn <- function(par) {
    # CASE A: If user provided a numeric tvec, use it directly
    # The only issue with this branch is that we hope the resulting tape is used only once.
    # It's addition here is to help with standard error calculations; when MLEs.tvec has been obtained.
    if (!is.null(user_tvec)) {
      # We still call `tvec_hat_for_ad()` to ensure the AD engine 
      # sees tvec as part of the expression.
      final_tvec <- tvec_hat_for_ad(
        theta         = par,
        tvec          = user_tvec,
        observations  = observed.data,
        K1_fn         = cgf$K1,
        K2_solve_fn   = K2_solve_fn
      )
      return(chosen_negll(final_tvec, par))
    }
    
    # CASE B: If cgf has an analytic formula for tvec_hat, use it
    if (cgf$has_analytic_tvec_hat()) {
      tvec_vals <- cgf$analytic_tvec_hat(observed.data, par)
      return(chosen_negll(tvec_vals, par))
    }
    
    # CASE C: numeric solve each time
    #   We do a numeric solve with the *numeric* version of par: RTMB:::getValues(par)
    #   That means AD won't see the iterative steps, but we do get a different tvec
    #   for each par. Then we pass it to tvec_hat_for_ad for AD tracking.
    ##### One potential issue is that getValues(par) is an internal RTMB function => it may change
    tvec_vals <- saddlepoint.solve(  
      theta = RTMB:::getValues(par),  # numeric version of par
      y     = observed.data,
      cgf   = cgf,
      ...
    )
    # tvec_vals is numeric
    final_tvec <- tvec_hat_for_ad(
      theta         = par,
      tvec          = tvec_vals,
      observations  = observed.data,
      K1_fn         = cgf$K1,
      K2_solve_fn   = K2_solve_fn
    )
    chosen_negll(final_tvec, par)
  }
  tape_obj <- MakeTape(f = local_fn, x = param_vec)
  

  if (hessian) jacfun_obj <- tape_obj$jacfun()
  
  function(par) {
    list(
      vals     = tape_obj(par),
      gradient = if(gradient) as.vector(tape_obj$jacobian(par)) else NULL,
      hessian  = if(hessian) jacfun_obj$jacobian(par) else NULL
    )
  }
  
}







#' Compute a saddlepoint negative log-likelihood
#'
#' @description
#' Computes either "standard" or "zeroth-order" saddlepoint approximation to the negative log-likelihood.
#'
#' @param parameter_vector Numeric vector of parameters.
#' @param observed.data Numeric vector of observed data.
#' @param cgf A `CGF` object.
#' @param tvec.hat (Optional) numeric vector. If provided, it overrides other ways
#'        of finding `tvec`.
#' @param gradient Logical. If `TRUE`, the gradient is computed.
#' @param hessian Logical. If `TRUE`, the Hessian is computed.
#' @param method Character string. One of `"standard"` or `"zeroth"`. 
#' @param ... Additional arguments passed to `saddlepoint.solve()` if `tvec.hat` is not provided.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{vals}{Numeric scalar. The negative log-likelihood at `parameter_vector`.}
#'   \item{gradient}{Numeric vector. The gradient at `parameter_vector`, if `gradient=TRUE`.}
#'   \item{hessian}{The Hessian at `parameter_vector`, if `hessian=TRUE`.}
#' }
#' 
#' @seealso \code{\link{saddlepoint.solve}}
#'
#' @export
compute.spa.negll <- function(parameter_vector,
                              observed.data,
                              cgf,
                              tvec.hat = NULL,
                              gradient = FALSE,
                              hessian  = FALSE,
                              method   = "standard",
                              ...) {
  if (!inherits(cgf, "CGF")) stop("`cgf` must be an object of class CGF")
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data))    stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be a logical.")
  if (!is.logical(hessian) || length(hessian) != 1) stop("`hessian` must be a logical.")
  if (!is.character(method) || length(method) != 1) stop("`method` must be a character string.")
  
  
  
  
  
  # If no gradient/hessian needed => direct numeric eval
  if (!gradient && !hessian) {
    final_tvec <- get_tvec_hat_vals(
      parameter_vector = parameter_vector,
      observed.data    = observed.data,
      cgf              = cgf,
      user_tvec        = tvec.hat,
      ...
    )
    chosen_negll <- choose_negll_fun(method = method, cgf = cgf)
    val <- chosen_negll(final_tvec, parameter_vector)[1] ##### some attribute has not been stripped off
    return(list(vals = val, gradient = NULL, hessian = NULL))
  }


  # Otherwise, get the taped function
  fn <- get_spa_taped_fun(
          param_vec     = parameter_vector,
          observed.data = observed.data,
          cgf           = cgf,
          method        = method,
          user_tvec     = tvec.hat,
          gradient      = gradient,
          hessian       = hessian,
          ...
  )
  fn(parameter_vector)
}



