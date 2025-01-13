# R/compute-funcT.R
# Objects: compute.funcT


# This script sets up the function(s) that compute the correction terms
# for the zeroth/first-order saddlepoint approximation to the log-likelihood.
# 
# The correction for the zeroth order is: -0.5*log(det(K2_val))
# The first-order's correction uses the private method func_T from the CGF object.



#' Compute the correction term for the saddlepoint approximation to the log-likelihood
#' 
#' @description
#' This function calculates the correction term for the saddlepoint approximation 
#' to the log-likelihood. By default, it computes that of the **first-order** approximation,
#' If \code{zeroth.order=TRUE}, it instead computes the correction term for the **zeroth-order** approximation.
#' 
#' @param theta Numeric vector of parameters.
#' @param observed.data Numeric vector of observed data.
#' @param cgf A CGF object.
#' @param tvec.hat (Optional) Numeric vector. If supplied, \code{tvec} is taken directly 
#'   as this vector (the saddlepoint \eqn{\hat{t}}). Otherwise, we solve for it or try 
#'   \code{analytic_tvec_hat}.
#' @param gradient Logical. If `TRUE`, the gradient wrt \code{theta} is returned.
#' @param hessian Logical. If `TRUE`, the Hessian wrt \code{theta} is returned. 
#' @param zeroth.order Logical. Default is `FALSE`. 
#' @param ... Optional arguments to be passed to \code{\link{saddlepoint.solve}}.
#' 
#' @seealso \code{\link{find.saddlepoint.MLE}}
#' 
#' @return A named list with elements:
#' \describe{
#'   \item{val}{Numeric vector of the correction term.}
#'   \item{gradient}{The gradient at `theta`, if `gradient=TRUE`.}
#'   \item{hessian}{The Hessian at `theta`, if `hessian=TRUE`.}
#' }
#' 
#' @export
compute.funcT <- function(theta, 
                          observed.data, 
                          cgf, 
                          tvec.hat = NULL, 
                          gradient = FALSE, 
                          hessian = FALSE, 
                          zeroth.order = FALSE, 
                          ...) {
  # Basic chcks
  if (!is.numeric(theta)) stop("'theta' is not defined for type ", class(theta))
  if (!is.numeric(observed.data)) stop("'observed.data' is not defined for type ", class(observed.data))
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("'gradient' must be a logical.")
  if (!is.logical(hessian) || length(hessian) != 1) stop("'hessian' must be a logical.")
  
  if (!is.logical(zeroth.order) || length(zeroth.order) != 1) stop("'zeroth.order' must be a single logical.")
  
  
  # Access the first-order correction method (func_T), if needed
  funcT_firstorder <- cgf$.get_private_method("func_T") # arguments: (tvec, theta)
  
  # helper for the zeroth-order correction
  zeroth_order_correction <- function(tvec, par_vec) {
    K2_val <- cgf$K2(tvec, par_vec)
    -0.5 * determinant(K2_val, logarithm = TRUE)$modulus
  }
  
  # Choose which correction function to use
  correction_func <- if (zeroth.order) {
    zeroth_order_correction
  } else {
    funcT_firstorder
  }
  
  # MakeTape requires a function of param only, 
  # so we define a second-level wrapper if we do AD
  # but that depends on how we do tvec_hat solving, etc.
  
  
  # Generate the dfdu_solve_fn function required by the C++ tvec_hat()
  # This function enables proper taping by capturing the dependency of theta through tvec
  dfdu_solve_fn <- create_tvec_hat_dfdu_solve_fn(cgf)
  
  # Check if gradient and/or Hessian is required
  taping_required <- (gradient || hessian)
  if (hessian && !gradient) gradient <- TRUE
  
  if (!is.null(tvec.hat)) {
    if (length(tvec.hat) != length(observed.data)) stop("`tvec.hat` must match the length of `observed.data`.")
    
    if (!taping_required) { return(list(vals = correction_func(tvec.hat, theta))) }
    
    func_T_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- tvec_hat(theta          = par_vec,
                                  tvec           = tvec.hat,
                                  observations   = observed.data,
                                  K1_fn          = cgf$K1,
                                  dfdu_solve_fn  = dfdu_solve_fn
      )
      correction_func(tvec.hat.for.ad, par_vec)
    }
    tape_obj <- MakeTape(f = func_T_wrapper, x = theta)
    return(get_tape_res(tape_obj, theta, gradient, hessian))
  }
  
  tvec.hat.analytic_res <- cgf$analytic_tvec_hat(observed.data, theta)
  if (!is.null(tvec.hat.analytic_res)) {
    
    if (!taping_required) { return( list(vals = correction_func(tvec.hat.analytic_res, theta) ) ) }
    
    func_T_wrapper <- function(par_vec) {
      tvec.hat.for.ad <- cgf$analytic_tvec_hat(observed.data, par_vec)
      correction_func(tvec.hat.for.ad, par_vec)
    }
    tape_obj <- MakeTape(f = func_T_wrapper, x = theta)
    return(get_tape_res(tape_obj, theta, gradient, hessian))
  }
  
  
  saddlepoint_res <- saddlepoint.solve(
    theta = theta,
    y     = observed.data,
    cgf   = cgf,
    ...
  )
  
  if (!taping_required) { return(list(vals = correction_func(saddlepoint_res, theta) )) }
  
  func_T_wrapper <- function(par_vec) {
    tvec.hat.for.ad <- tvec_hat(theta          = par_vec,
                                tvec           = saddlepoint_res,
                                observations   = observed.data,
                                K1_fn          = cgf$K1,
                                dfdu_solve_fn  = dfdu_solve_fn
    )
    correction_func(tvec.hat.for.ad, par_vec)
  }
  tape_obj <- MakeTape(f = func_T_wrapper, x = theta)
  get_tape_res(tape_obj, theta, gradient, hessian)
  
  
}