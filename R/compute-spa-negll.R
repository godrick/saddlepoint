# R/compute-spa-negll.R
# Objects: compute.spa.negll, compute.zeroth.order.spa.negll, get_tape_res 




#-------------------------------------------------------------------
# Helper for returning results
#-------------------------------------------------------------------
#' @noRd
get_tape_res <- function(tape_obj, param_vec, gradient = FALSE, hessian = FALSE) {
  
  out_list <- list(
    vals = tape_obj(param_vec),
    gradient = NULL,
    hessian = NULL
  )
  
  if (gradient) { out_list$gradient <- as.vector(tape_obj$jacobian(param_vec)) }
  
  if (hessian) {
    jacfun_obj <- tape_obj$jacfun()
    out_list$hessian <- jacfun_obj$jacobian(param_vec)
  }
  out_list
}


#-------------------------------------------------------------------
# Main Functions
#-------------------------------------------------------------------

#' Compute saddlepoint approximation-based negative log-likelihood
#'
#' @description
#' Computes the saddlepoint-based negative log-likelihood for given parameters and data using one of the these:
#' - **User-Provided `tvec.hat`:** If `tvec.hat` is supplied, it is used directly.
#' - **Analytic `tvec.hat`:** If the `CGF` object provides an analytic solution via `analytic_tvec_hat()`, it is utilized.
#' - **Saddlepoint Solver:** Otherwise, a numeric solution is obtained using `saddlepoint.solve()`.
#'
#' Internally, this function creates a tape (via RTMB's `MakeTape`) for the negative log-likelihood with respect to the parameter vector, 
#' and returns the evaluated negative log-likelihood, the gradient, and optionally the Hessian.
#'
#' @param parameter_vector Numeric vector of parameters.
#' @param observed.data    Numeric vector of observed data.
#' @param cgf          A `CGF` object. 
#' @param tvec.hat    (Optional) Numeric vector. If supplied, `tvec` is taken directly as this vector.
#' @param gradient    Logical. If `TRUE`, the gradient is computed; if `FALSE`, it is not returned. Defaults to `FALSE`.
#' @param hessian     Logical. If `TRUE`, the Hessian is computed; if `FALSE`, a `NULL` is returned in its place.
#' @param ...         Further arguments passed to `saddlepoint.solve()` (if needed).
#' 
#' @seealso \code{\link{saddlepoint.solve}}, \code{\link{compute.zeroth.order.spa.negll}}
#'
#' @return A named list with elements:
#' \describe{
#'   \item{vals}{Numeric scalar. The negative log-likelihood at `parameter_vector`.}
#'   \item{gradient}{Numeric vector. The gradient at `parameter_vector`, if `gradient=TRUE`.}
#'   \item{hessian}{The Hessian at `parameter_vector`, if `hessian=TRUE`.}
#' }
#' 
#' @export
compute.spa.negll <- function(parameter_vector,
                              observed.data,
                              cgf,
                              tvec.hat = NULL,
                              gradient = FALSE,
                              hessian = FALSE,
                              ...) {
  
  neg_ll <- cgf$.get_private_method("neg_ll")
  
  # Validate arguments
  if (!inherits(cgf, "CGF")) {
    stop("`cgf` must be an object of class CGF")
  }
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data))    stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be a logical.")
  if (!is.logical(hessian) || length(hessian) != 1) stop("`hessian` must be a logical.")
  
  
  
  # Generate the dfdu_solve_fn function required by the C++ tvec_hat()
  # This function enables proper taping by capturing the dependency of theta through tvec
  dfdu_solve_fn <- create_tvec_hat_dfdu_solve_fn(cgf)
  
  
  # Check if gradient and/or Hessian is required
  taping_required <- gradient || hessian
  if (hessian && !gradient) gradient <- TRUE
  
  
  #-------------------------------------------------------------------
  # If user already provided tvec => skip any solving
  # We might use this for standard error functions where optimal tvec is known.
  #-------------------------------------------------------------------
  if (!is.null(tvec.hat)) {
    if (length(tvec.hat) != length(observed.data)) stop("`tvec.hat` must match the length of `observed.data`.")
    
    if (!taping_required) { return( list(vals = neg_ll(tvec.hat, parameter_vector)) ) }
    
    neg_ll_spa_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- tvec_hat(theta          = par_vec,
                                  tvec           = tvec.hat,
                                  observations   = observed.data,
                                  K1_fn          = cgf$K1,
                                  # K2_fn        = cgf$K2
                                  dfdu_solve_fn  = dfdu_solve_fn
                                  )
      neg_ll(tvec.hat.for.ad, par_vec)
    }
    tape_obj <- MakeTape(f = neg_ll_spa_wrapper, x = parameter_vector)
    return(get_tape_res(tape_obj, parameter_vector, gradient, hessian))
  }
  
  
  #-------------------------------------------------------------------
  # If `cgf$analytic_tvec_hat(...)` returns a valid numeric vector
  #-------------------------------------------------------------------
  tvec.hat.analytic_res <- cgf$analytic_tvec_hat(observed.data, parameter_vector)
  if (!is.null(tvec.hat.analytic_res)) {
    
    if (!taping_required) { return( list(vals = neg_ll(tvec.hat.analytic_res, parameter_vector)) ) }
    
    neg_ll_spa_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- cgf$analytic_tvec_hat(observed.data, par_vec)
      neg_ll(tvec.hat.for.ad, par_vec)
    }
    tape_obj <- MakeTape(f = neg_ll_spa_wrapper, x = parameter_vector)
    return(get_tape_res(tape_obj, parameter_vector, gradient, hessian))
  }
  
  #-------------------------------------------------------------------
  # Otherwise, solve via saddlepoint.solve(...)
  #-------------------------------------------------------------------
  saddlepoint_res <- saddlepoint.solve(
    theta = parameter_vector,
    y     = observed.data,
    cgf   = cgf,
    ...
  )
  
  if (!taping_required) { return(list(vals = neg_ll(saddlepoint_res, parameter_vector))) }
  
  neg_ll_spa_wrapper <- function(par_vec) {
    tvec.hat.for.ad <- tvec_hat(theta          = par_vec,
                                tvec           = saddlepoint_res,
                                observations   = observed.data,
                                K1_fn          = cgf$K1,
                                dfdu_solve_fn  = dfdu_solve_fn
                                )
    neg_ll(tvec.hat.for.ad, par_vec)
  }
  tape_obj <- MakeTape(f = neg_ll_spa_wrapper, x = parameter_vector)
  get_tape_res(tape_obj, parameter_vector, gradient, hessian)
}
















#' Compute Zeroth-order saddlepoint approximation-based negative log-likelihood
#'
#' @description
#' Computes the zeroth-order saddlepoint approximation-based negative log-likelihood for given parameters and data.
#' 
#'
#' The zeroth-order approximation is defined as the negative of the tilting exponent.
#'
#'
#' @param parameter_vector Numeric vector of parameters.
#' @param observed.data    Numeric vector of observed data.
#' @param cgf              A `CGF` object. Must provide methods `K1()`, `K2()`, and optionally `analytic_tvec_hat()`.
#' @param tvec.hat         (Optional) Numeric vector. If supplied, it is used directly as `tvec`. 
#' @param gradient         Logical. If `TRUE`, the gradient is computed; if `FALSE`, it is not returned. Defaults to `FALSE`.
#' @param hessian          Logical. If `TRUE`, the Hessian is computed; if `FALSE`, `hessian` is returned as `NULL`.
#' @param ...              Additional arguments passed to `saddlepoint.solve()` (if needed).
#' 
#' @seealso \code{\link{saddlepoint.solve}}, \code{\link{compute.spa.negll}}
#'
#' @return A named list with elements:
#' \describe{
#'   \item{vals}{Numeric scalar. The zeroth-order saddlepoint negative log-likelihood at `parameter_vector`.}
#'   \item{gradient}{Numeric vector. The gradient at `parameter_vector`, if `gradient=TRUE`.}
#'   \item{hessian}{The Hessian at `parameter_vector`, if `hessian=TRUE`.}
#' }
#' 
#'
#' @export
compute.zeroth.order.spa.negll <- function(parameter_vector,
                                           observed.data,
                                           cgf,
                                           tvec.hat = NULL,
                                           gradient = FALSE,
                                           hessian = FALSE,
                                           ...) {
  
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  
  
  if (!inherits(cgf, "CGF")) stop("`cgf` must be an object of class CGF")
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data))    stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be a logical.")
  if (!is.logical(hessian) || length(hessian) != 1) stop("`hessian` must be a logical.")
  
  # Generate the dfdu_solve_fn function required by the C++ tvec_hat()
  dfdu_solve_fn <- create_tvec_hat_dfdu_solve_fn(cgf)
  
  # Check if gradient and/or Hessian is required
  taping_required <- gradient || hessian
  if (hessian && !gradient) gradient <- TRUE
  
  #-------------------------------------------------------------------
  # If user provided tvec.hat, use it directly
  #-------------------------------------------------------------------
  if (!is.null(tvec.hat)) {
    
    if (length(tvec.hat) != length(observed.data)) stop("`tvec.hat` must match the length of `observed.data`.")
    if (!taping_required) { return( list(vals = -tilting_exponent(tvec.hat, parameter_vector)) ) }
    
    # Create the wrapper function for negative log-likelihood
    zeroth.neg_ll_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- tvec_hat(tvec         = tvec.hat,
                                  theta        = par_vec,
                                  observations = observed.data,
                                  K1_fn        = cgf$K1,
                                  dfdu_solve_fn  = dfdu_solve_fn
                                  )
      -tilting_exponent(tvec.hat.for.ad, par_vec)
    }
    
    tape_obj <- MakeTape(f = zeroth.neg_ll_wrapper, x = parameter_vector)
    return(get_tape_res(tape_obj, parameter_vector, gradient,  hessian))
  }
  
  #-------------------------------------------------------------------
  # If `cgf$analytic_tvec_hat` provides a valid tvec.hat
  #-------------------------------------------------------------------
  tvec.hat.analytic_res <- cgf$analytic_tvec_hat(observed.data, parameter_vector)
  if (!is.null(tvec.hat.analytic_res)) {
    
    if (!taping_required) { return( list(vals = -tilting_exponent(tvec.hat.analytic_res, parameter_vector)) ) }
    
    zeroth.neg_ll_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- cgf$analytic_tvec_hat(observed.data, par_vec)
      -tilting_exponent(tvec.hat.for.ad, par_vec)
    }
    
    tape_obj <- MakeTape(f = zeroth.neg_ll_wrapper, x = parameter_vector)
    return(get_tape_res(tape_obj, parameter_vector, gradient, hessian))
  }
  
  #-------------------------------------------------------------------
  # Otherwise, solve using saddlepoint.solve()
  #-------------------------------------------------------------------
  saddlepoint_res <- saddlepoint.solve(
    theta = parameter_vector,
    y     = observed.data,
    cgf   = cgf,
    ...
  )
  
  if (!taping_required) { return(list(vals = -tilting_exponent(saddlepoint_res, parameter_vector)) ) }
  
  zeroth.neg_ll_wrapper <- function(par_vec) {
    tvec.hat.for.ad <- tvec_hat(tvec         = saddlepoint_res,
                                theta        = par_vec,
                                observations = observed.data,
                                K1_fn        = cgf$K1,
                                dfdu_solve_fn  = dfdu_solve_fn
                                )
    -tilting_exponent(tvec.hat.for.ad, par_vec)
  }
  
  tape_obj <- MakeTape(f = zeroth.neg_ll_wrapper, x = parameter_vector)
  get_tape_res(tape_obj, parameter_vector, gradient, hessian)
}





