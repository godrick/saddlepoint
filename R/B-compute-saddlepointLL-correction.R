# R/B-compute-saddlepointLL-correction.R
# Objects: compute.saddlepointLL.correction


# This script sets up the function(s) that compute the correction terms
# for the zeroth/first-order saddlepoint approximation to the log-likelihood.
# 
# The correction for the zeroth order is: -0.5*log(det(K2_val))
# The first-order's correction uses the private method func_T from the CGF object.

#### Helper functions used can be found in the file: B-compute-spa-negll.R








#' Compute the correction term for the saddlepoint approximation to the log-likelihood
#' 
#' @description
#' This function calculates the correction term for the saddlepoint approximation 
#' to the log-likelihood. By default, it computes that of the **first-order** approximation,
#' If \code{spa_method="zeroth"}, it instead computes the correction term for the **zeroth-order** approximation.
#' 
#' @param parameter_vector Numeric vector of parameters.
#' @param observed.data Numeric vector of observed data.
#' @param cgf A `CGF` object.
#' @param tvec.hat (Optional) Numeric vector. If supplied, \code{tvec} is taken directly 
#'   as this vector (the saddlepoint \eqn{\hat{t}}). Otherwise, we solve for it numerically if \code{cgf$has_analytic_tvec_hat()} is `FALSE`.
#' @param gradient Logical. If `TRUE`, the gradient wrt \code{parameter_vector} is returned.
#' @param hessian Logical. If `TRUE`, the Hessian wrt \code{parameter_vector} is returned. 
#' @param spa_method Character string. One of `"standard"` or `"zeroth"`. 
#' @param ... Optional arguments to be passed to \code{\link{saddlepoint.solve}}.
#' 
#' @seealso \code{\link{find.saddlepoint.MLE}}, \code{\link{compute.spa.negll}}
#' 
#' @return A named list with elements:
#' \describe{
#'   \item{val}{Numeric scalar of the correction term.}
#'   \item{gradient}{The gradient at `theta`, if `gradient=TRUE`.}
#'   \item{hessian}{The Hessian at `theta`, if `hessian=TRUE`.}
#' }
#' 
#' @export
compute.saddlepointLL.correction <- function(parameter_vector,
                                             observed.data,
                                             cgf,
                                             tvec.hat = NULL,
                                             gradient = FALSE,
                                             hessian  = FALSE,
                                             spa_method   = "standard",
                                             ...) {
  if (!inherits(cgf, "CGF")) stop("`cgf` must be an object of class CGF")
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data)) stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be a single logical.")
  if (!is.logical(hessian) || length(hessian) != 1) stop("`hessian` must be a single logical.")
  if (!is.character(spa_method) || length(spa_method) != 1) stop("`spa_method` must be a single character string.")
  
  if (spa_method == "standard") {
    spa_method <- "correction_standard"
  } else if (spa_method == "zeroth") {
    spa_method <- "correction_zeroth"
  }
  
  if (!gradient && !hessian) {
    tvec_hat_vals <- get_nonAD_tvec_hat_vals(
      parameter_vector = parameter_vector,
      observed.data    = observed.data,
      cgf              = cgf,
      user_tvec        = tvec.hat,
      ...
    )
    chosen_negll <- choose_spa_function(spa_method = spa_method, cgf = cgf)
    val <- chosen_negll(tvec_hat_vals, parameter_vector)[1] ##### here, [1] strips off an attribute.
    return(list(vals = val, gradient = NULL, hessian = NULL))
  }
  
  # Otherwise, build an AD tape and return its evaluation.
  taped_fun <- create_spa_taped_fun(
    param_vec     = parameter_vector,
    observed.data = observed.data,
    cgf           = cgf,
    spa_method    = spa_method,
    user_tvec     = tvec.hat,
    gradient      = gradient,
    hessian       = hessian,
    ... #### addtional arguments are passed to saddlepoint.solve() are still not being passed to cpp
  )
  taped_fun(parameter_vector)
  
  
}



