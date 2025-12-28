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
    K_elem = function(tvec, pm) {
      -log(pm[,1] - tvec) + log(pm[,1])
    },
    K1_elem     = function(tvec, pm) {
      # print(class(pm))
      # print(class(tvec))
      1/(pm[,1] - tvec)
    },
    # K2_elem     = function(tvec, pm)  1 / (pm[,1] - tvec)^2,
    K2_elem     = function(tvec, pm)  {
      denom <- pm[,1] - tvec
      denom2 <- denom * denom
      1 / denom2
    },
    # K3_elem     = function(tvec, pm)  2 / (pm[,1] - tvec)^3,
    K3_elem     = function(tvec, pm)  {
      denom <- pm[,1] - tvec
      denom2 <- denom * denom
      denom3 <- denom2 * denom
      2 / denom3
    },
    # K4_elem     = function(tvec, pm)  6 / (pm[,1] - tvec)^4,
    K4_elem     = function(tvec, pm)  {
      denom <- pm[,1] - tvec
      denom2 <- denom * denom
      denom4 <- denom2 * denom2
      6 / denom4
    },
    t_hat_elem  = function(x, pm) pm[,1] - 1/x,
    ineq_elem = function(tvec, pm) {
      # print(class(tvec))
      # print(class(pm))
      tvec - pm[,1]
    },
    split_param_to_mat = function(param) {
      cbind(param)
      # matrix(param, ncol = 1L)
    },
    simulate_func = function(iidReps, parameter_vector, ...) {
      rate <- as.numeric(parameter_vector)

      if (any(!is.finite(rate)) || any(rate <= 0)) stop("ExponentialCGF$rsim: 'rate' must be finite and > 0.")


      d <- length(rate)
      out <- stats::rexp(
        n    = d * iidReps,
        rate = rep.int(rate, times = iidReps)
      )
      matrix(out, nrow = d, ncol = iidReps)
    },

    iidReps = iidReps,
    op_name = op_name,
    ...
  )
}



# .exponential_base_cgf <- function(iidReps, op_name, ...) {
#   .make_univariate_model_cgf_matrix(
#     K_elem      = function(tvec, pm) -log1p(-tvec / pm[,1]),  # stable form
#     K1_elem     = function(tvec, pm)  1 / (pm[,1] - tvec),
#     K2_elem     = function(tvec, pm)  1 / (pm[,1] - tvec)^2,
#     K3_elem     = function(tvec, pm)  2 / (pm[,1] - tvec)^3,
#     K4_elem     = function(tvec, pm)  6 / (pm[,1] - tvec)^4,
#     t_hat_elem  = function(x, pm) pm[,1] - 1 / x,
#     split_param_to_mat = function(param) matrix(param, ncol = 1L),
#     iidReps = iidReps,
#     op_name = op_name,
#     ineq_elem = function(tvec, pm) tvec - pm[,1],
#     ...
#   )
# }





#' Exponential CGF Object
#'
#' A ready-to-use CGF object for the Exponential distribution.
#' This object is vectorized for i.i.d. replicates of the same rate.
#'
#' @format An object of class \code{CGF} (R6), with usual methods \code{K},
#'   \code{K1}, \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' ExponentialCGF$K1(0, 2)  # (expected value for Exp(rate=2) is 1/2)
#'
#' @export
ExponentialCGF <- .exponential_base_cgf(iidReps = "any", op_name = "ExponentialCGF")

# ----------------------------------------------------------------------------
#   A Parametric Model CGF for the Exponential distribution.
# ----------------------------------------------------------------------------






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
  base_cgf <- .exponential_base_cgf(iidReps = iidReps, op_name = "ExponentialModelCGF", ...)
  adaptCGF(cgf = base_cgf, adaptor = rate_fn)
}











