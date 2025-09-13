# R/ExponentialCGF.R
# Objects:









##############################################################################
## Exponential CGF Setup
##############################################################################
##
## Distribution: Exponential(rate)
## CGF: K(t) = -log(rate - t) + log(rate),  valid for t < rate.
##
##############################################################################

# ----------------------------------------------------------------------------
# First, a "ready-to-use" CGF object for a single rate parameter.
# ----------------------------------------------------------------------------





## Internal factory – build Exponential CGF via univariate utilities
.exponential_base_cgf <- function(iidReps, op_name, ...) {
  .make_univariate_model_cgf_matrix(
    K_elem      = function(tvec, pm) -log(pm[,1] - tvec) + log(pm[,1]),
    K1_elem     = function(tvec, pm)  1 / (pm[,1] - tvec),
    K2_elem     = function(tvec, pm)  1 / (pm[,1] - tvec)^2,
    K3_elem     = function(tvec, pm)  2 / (pm[,1] - tvec)^3,
    K4_elem     = function(tvec, pm)  6 / (pm[,1] - tvec)^4,
    t_hat_elem  = function(x, pm) pm[,1] - 1 / x,
    split_param_to_mat = function(param) matrix(param, ncol = 1L),
    iidReps = iidReps,
    op_name = op_name,
    ineq_elem = function(tvec, pm) tvec - pm[,1],
    ...
  )
}

#' Exponential CGF Object
#'
#' A ready-to-use CGF object for the Exponential distribution.
# This object is vectorized for i.i.d. replicates of the same rate.
#' 
#' @format An object of class \code{CGF} (R6), with usual methods \code{K},
#'   \code{K1}, \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' ExponentialCGF$K1(0, 2)  # (expected value for Exp(rate=2) is 1/2)
#'
#' @export
ExponentialCGF <- .exponential_base_cgf(iidReps = "any", op_name = "ExponentialCGF")
  K_vectorized_func = function(tvec, lambda) {
    -log(lambda[1] - tvec) + log(lambda[1])
  },
  K1_vectorized_func = function(tvec, lambda) { 1 / (lambda[1] - tvec) }, 
  K2_vectorized_func = function(tvec, lambda) { 1 / (lambda[1] - tvec)^2},
  K3_vectorized_func = function(tvec, lambda) { 2 / (lambda[1] - tvec)^3},
  K4_vectorized_func = function(tvec, lambda) { 6 / (lambda[1] - tvec)^4},
  ineq_constraint_func = function(tvec, lambda) { 
    # Enforce tvec < lambda => tvec - lambda < 0
    tvec - lambda[1] 
  },
  analytic_tvec_hat_func = function(x, lambda) {
    # 1/(lambda - t)= x => t= lambda - 1/x
    lambda[1] - 1 / x
  },
  op_name = "ExponentialCGF"
)

# ----------------------------------------------------------------------------
#   A Parametric Model CGF for the Exponential distribution.
# ----------------------------------------------------------------------------


#' @noRd
validateExponentialLengths <- function(vec, lambda, iidReps) {
  len_vec <- length(vec)
  len_lambda <- length(lambda)
  if (!is.null(iidReps)) {
    expected <- len_lambda * iidReps
    if (len_vec != expected) {
      stop(sprintf("Length mismatch: input vector has length %d; expected %d (lambda length %d times iidReps %d).",
                   len_vec, expected, len_lambda, iidReps))
    }
  } else if (len_vec %% len_lambda != 0) {
    stop(sprintf("Length mismatch: input vector length %d is not a multiple of lambda length %d.",
                 len_vec, len_lambda))
  }
}

#' @noRd
.ExponentialModelCGF_internal <- function(iidReps, ...) {
  .exponential_base_cgf(iidReps = iidReps, op_name = "ExponentialModelCGF", ...)
}

#' Create a Parametric Exponential CGF Object
#'
#' @description
#' Constructs a CGF object for the Exponential distribution where the rate parameter
#' \eqn{\lambda(\theta)} is derived from a user-supplied function (or adaptor).
#' This supports both i.i.d. and non-identical usage, depending on whether \eqn{\lambda(\theta)}
#' returns one or multiple values, and depending on the \code{iidReps} setting.
#'
#' @param rate A function (or adaptor) that accepts a parameter vector \code{theta} and returns
#'   the rate parameter \eqn{\lambda} (a positive numeric value or vector).
#' @param iidReps Either \code{"any"} or a positive integer specifying
#'   the number of i.i.d. blocks are expected. Each block correspond to one copy of the Exponential variables (or multiple if \eqn{\lambda(\theta)} is a vector).
#' @param ... Additional arguments passed to the underlying CGF creation function.
#'
#' @return A `CGF` object.
#'
#' @examples
#' rate_func <- function(theta) theta[1]  # For example, theta -> 2 gives rate = 2
#' expo_model_cgf <- ExponentialModelCGF(rate = rate_func, iidReps = 1)
#' # Evaluate the first derivative at t = 0 for rate 2:
#' expo_model_cgf$K1(0, c(2))
#'
#' @export
ExponentialModelCGF <- function(rate, iidReps = "any", ...) {
  .check_iidReps(iidReps)
  rate_fn <- validate_function_or_adaptor(rate)
  base_cgf <- .ExponentialModelCGF_internal(iidReps, ...)
  adaptCGF(cgf = base_cgf, adaptor = rate_fn)
}












