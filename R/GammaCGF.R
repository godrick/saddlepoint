# R/GammaCGF.R
# Objects: GammaCGF, GammaModelCGF

#' Gamma CGF Object
#'
#' A ready-to-use CGF object for the Gamma distribution with shape \eqn{\alpha}
#' and rate \eqn{\beta}. The \code{parameter_vector} used when calling methods like `K(tvec, parameter_vector)`
#' shoulde be a numeric vector \eqn{c(\alpha, \beta)}.
#'
#' @details
#' **CGF**: For a Gamma random variable \eqn{X} with shape \eqn{\alpha} and rate
#' \eqn{\beta}, the cumulant generating function is:
#' \deqn{K(t;\alpha, \beta) = -\alpha \,\log \bigl(1 - t/\beta\bigr), \quad t < \beta.}
#'
#' **Parameter Vector**: The \code{parameter_vector} is assumed to have the form 
#' \eqn{(\alpha, \beta)}. You must ensure 
#' that \eqn{t < \beta} for all evaluations. This object enforces that constraint 
#' as an inequality \eqn{t - \beta < 0}.
#'
#' @format An object of class \code{CGF} (an R6 class), with the usual methods:
#' \code{K, K1, K2, K3operator, K4operator}, etc.
#'
#' @examples
#' # Evaluate K at t=0.5 for shape=2, rate=2 (thus t<2).
#' # param = c(2, 2)
#' GammaCGF$K(0.5, c(2,2))
#'
#' @export
GammaCGF <- createCGF_fromVectorisedFunctions(
  
  K_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K(t) = -alpha * log1p(-t/beta)
    # using log1p(x) for x = -(tvec/beta)
    -alpha * log1p(-tvec/beta)
  },
  
  ########################
  ## K1(t)
  ########################
  K1_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K1(t) = alpha / (beta - t)
    alpha / (beta - tvec)
  },
  
  
  K2_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K2(t) = alpha / (beta - t)^2
    alpha / (beta - tvec)^2
  },
  
  
  K3_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    2 * alpha / (beta - tvec)^3
  },
  
  K4_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    6 * alpha / (beta - tvec)^4
  },
  
  
  # We want tvec < beta => tvec - beta < 0
  # We'll just return (tvec - beta), which must be negative
  ineq_constraint = function(tvec, param) {
    beta <- param[2]
    tvec - beta
  },
  
  analytic_tvec_hat_func = function(y, param)  { param[2] - (param[1] / y) },
  
  param_adaptor = function(x) x[1:2],
  op_name       = "GammaCGF"
)


#' Create a Parametric Gamma CGF Object
#'
#' @description
#' Creates a CGF object for the Gamma distribution with shape \eqn{\alpha(\theta)} and 
#' rate \eqn{\beta(\theta)} defined by user-provided parameter functions. 
#' The user supplies a function that maps their full parameter vector 
#' \code{theta} to \eqn{(\alpha, \beta)}.
#'
#' @param shape A function that accepts a single parameter vector \code{theta} and returns the shape parameter.
#' @param rate A function that accepts a single parameter vector \code{theta} and returns the rate parameter.
#' @param iidReps A positive integer specifying the number of i.i.d. replicates. Defaults to \code{1}.
#' @param ... Additional arguments passed to the underlying CGF creation function. 
#'   This can include optional overrides to the base CGF class.
#'
#'
#' @return A CGF object.
#' @export
GammaModelCGF <- function(shape,
                          rate,
                          iidReps = 1,
                          ...) {
  # Check user-supplied function
  if (!is.function(shape)) stop("`shape` must be a function.")
  if (length(formals(shape)) != 1) stop("`shape` must be a function that accepts exactly one argument.")
  if (!is.function(rate)) stop("`rate` must be a function")
  if (length(formals(rate)) != 1) stop("`rate` must be a function that accepts exactly one argument.")
  
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) stop(sprintf("`tvec` must have length %d (got %d).", iidReps, length(tvec)))
  }
  
  
  
  # param_adaptor => shape_rate
  base_cgf <- createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      -alpha * log1p(-tvec / beta)
    },
    K1_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      alpha / (beta - tvec)
    },
    K2_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      alpha / (beta - tvec)^2
    },
    K3_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      2*alpha / (beta - tvec)^3
    },
    K4_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      6*alpha / (beta - tvec)^4
    },
    ineq_constraint_func = function(tvec, sr) {
      check_tvec(tvec)
      beta <- sr[2]
      tvec - beta
    },
    param_adaptor = function(x) c(shape(x), rate(x)),
    analytic_tvec_hat_func = function(y, sr)  {
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      alpha <- sr[1]
      beta  <- sr[2]
      beta - (alpha / y)
    },
    op_name = "GammaModelCGF",
    ...
  )
  
  base_cgf
}

#### TODO: GammaNonIdenticalModelCGF



