# R/B-compute-saddlepointLL-correction.R
# Objects: compute.saddlepointLL.correction


# Correction terms for the SPA log-likelihood:
#   - Zeroth-order  :  -0.5 * log det K2(t, theta)
#   - First-order T :  uses cgf private method func_T (already AD-safe in base)

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
#'   as this vector (the saddlepoint \eqn{\hat{t}}). Otherwise we compute it (analytic if
#'   available, else numeric solve).
#' @param gradient Logical. If `TRUE`, return gradient wrt \code{parameter_vector}.
#' @param hessian Logical. If `TRUE`, return Hessian wrt \code{parameter_vector}.
#' @param spa_method Character string. One of `"standard"` (first-order) or `"zeroth"`.
#' @param tvec_source One of "auto","user","analytic","newton","solver_atomic". See `compute.spa.negll`.
#' @param solver_fun Numeric solver used when \code{tvec_source="solver_atomic"} or in non-AD path. Default `saddlepoint.solve`.
#' @param newton_t_init Optional start for t in the Newton solver path.
#' @param ... Passed to \code{solver_fun} when used.
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
                                             tvec.hat     = NULL,
                                             gradient     = FALSE,
                                             hessian      = FALSE,
                                             spa_method   = "standard",
                                             tvec_source  = c("auto","user","analytic","newton","solver_atomic"),
                                             solver_fun   = saddlepoint.solve,
                                             newton_t_init = NULL,
                                             ...) {
  if (!inherits(cgf, "CGF")) stop("`cgf` must be an object of class CGF")
  if (!is.numeric(parameter_vector)) stop("`parameter_vector` must be numeric.")
  if (!is.numeric(observed.data))    stop("`observed.data` must be numeric.")
  if (!is.null(tvec.hat) && !is.numeric(tvec.hat)) stop("`tvec.hat` must be numeric.")
  if (!is.logical(gradient) || length(gradient) != 1) stop("`gradient` must be logical(1).")
  if (!is.logical(hessian)  || length(hessian)  != 1) stop("`hessian` must be logical(1).")


  if (spa_method == "standard") {
    spa_method <- "correction_standard"
  } else if (spa_method == "zeroth") {
    spa_method <- "correction_zeroth"
  }

  tvec_source <- match.arg(tvec_source)

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
    chosen <- choose_spa_function(spa_method = spa_method, cgf = cgf)
    val <- chosen(tvec_hat_vals, parameter_vector)[1]
    return(list(vals = val, gradient = NULL, hessian = NULL))
  }

  taped_fun <- create_spa_taped_fun(
    param_vec     = parameter_vector,
    observed.data = observed.data,
    cgf           = cgf,
    spa_method    = spa_method,  # one of "correction_*"
    tvec_source   = if (!is.null(tvec.hat) && tvec_source == "auto") "user" else tvec_source,
    user_tvec     = tvec.hat,
    gradient      = gradient,
    hessian       = hessian,
    solver_fun    = solver_fun,
    newton_t_init = newton_t_init,
    ...
  )
  taped_fun(parameter_vector)
}




