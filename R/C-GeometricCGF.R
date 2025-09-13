# R/C-GeometricCGF.R

##############################################################################
## Geometric CGF Setup
##############################################################################
##
## Distribution: Geometric(p), counting the number of failures before first success.
## PMF: P(X = k) = (1 - p)^k * p,  k = 0,1,2,...
## 
## CGF reference: 
##    K(t; p) = log( p ) - log( 1 - exp(t) + p * exp(t) ),
## valid for t < -log(1 - p). 
##
##############################################################################

# ----------------------------------------------------------------------------
# First: A "ready-to-use" CGF object for a single geometric parameter "p".
#        The user calls this with the scalar parameter_vector = p. By default, this allows i.i.d. usage. 
#        For non-identical usage, see the second approach below.
# ----------------------------------------------------------------------------


## Internal factory – build Geometric CGF via univariate utilities
.geometric_base_cgf <- function(iidReps, op_name, ...) {
  .make_univariate_model_cgf_matrix(
    K_elem      = function(tvec, pm) log(pm[,1]) - log(1 - exp(tvec) + pm[,1]*exp(tvec)),
    K1_elem     = function(tvec, pm) (exp(tvec) - pm[,1]*exp(tvec)) / (1 - exp(tvec) + pm[,1]*exp(tvec)),
    K2_elem     = function(tvec, pm) {
      tmp <- 1 - exp(tvec) + pm[,1]*exp(tvec)
      (exp(tvec) - pm[,1]*exp(tvec)) / tmp^2
    },
    K3_elem     = function(tvec, pm) {
      tmp <- 1 - exp(tvec) + pm[,1]*exp(tvec)
      (exp(tvec) - pm[,1]*exp(tvec)) * (1 + exp(tvec) - pm[,1]*exp(tvec)) / tmp^3
    },
    K4_elem     = function(tvec, pm) {
      tmp <- 1 - exp(tvec) + pm[,1]*exp(tvec)
      (exp(tvec) - pm[,1]*exp(tvec)) * (1 + exp(2*tvec) + 4*exp(tvec)
        - 2*pm[,1]*exp(2*tvec) - 4*pm[,1]*exp(tvec) + pm[,1]^2*exp(2*tvec)) / tmp^4
    },
    t_hat_elem  = function(x, pm) log(x) - log(1 + x - pm[,1] - pm[,1]*x),
    split_param_to_mat = function(param) matrix(param, ncol = 1L),
    iidReps = iidReps,
    op_name = op_name,
    ineq_elem = function(tvec, pm) (1 - pm[,1]) * exp(tvec) - 1,
    ...
  )
}

#' Geometric CGF Object
#'
#' A ready-to-use CGF object for a single-parameter Geometric distribution.
#' This corresponds to the count of failures before the first success.
#' By default, this object is vectorized for i.i.d. replicates of probability `prob`.
#'
#' @seealso \code{\link{GeometricModelCGF}}
#'
#' @format An object of class \code{CGF} (R6), with usual methods:
#' \code{K, K1, K2, K3operator, K4operator}, etc.
#'
#' @examples
#' # expected value of X~Geometric(prob = 0.3) via CGF
#' GeometricCGF$K1(0, 0.3)
#'
#' @export
GeometricCGF <- .geometric_base_cgf(iidReps = "any", op_name = "GeometricCGF")











# ----------------------------------------------------------------------------
# Next, A parametric model, for i.i.d. and non-identical usage.
# ----------------------------------------------------------------------------


.GeometricModelCGF_internal <- function(iidReps, ...){
  .geometric_base_cgf(iidReps = iidReps, op_name = "GeometricModelCGF", ...)
}
















#' Create a Parametric Geometric CGF Object
#'
#' @description
#' Creates a CGF object for the Geometric distribution with success probability
#' \eqn{prob(\theta)} given by a user-supplied function or adaptor.
#' The resulting object correspond to the random variable that counts the number of failures before achieving the first success.
#' It supports i.i.d. and non-identical contexts with optional length enforcement via `iidReps`.
#'
#' @param prob A function (or adaptor) that accepts a single parameter vector \code{theta}
#'   and returns the success probability \eqn{prob} (a scalar) or a vector of probabilities.
#' @param iidReps Either \code{"any"} (no forced dimension) or a positive integer specifying how many
#'   i.i.d. blocks are expected. Each block correspond to one copy of the geometric variables (or multiple if `prob` is a vector).
#' @param ... Additional arguments passed to the underlying CGF creation function
#'
#' @return A `CGF` object
#' @export
GeometricModelCGF <- function(prob, iidReps = "any", ...) {
  .check_iidReps(iidReps)
  p_fn <- validate_function_or_adaptor(prob)
  base_cgf <- .GeometricModelCGF_internal(iidReps, ...)
  adaptCGF(cgf = base_cgf, adaptor = p_fn)
}
