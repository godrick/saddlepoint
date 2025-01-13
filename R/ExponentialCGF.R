# R/ExponentialCGF.R

#' Exponential CGF Object
#'
#' A ready-to-use CGF object for the Exponential distribution.
#'
#' @details
#' **Parameter Vector**:
#' The parameter vector must be a single positive numeric value representing the rate parameter.
#'
#' **Cumulant Generating Function**:
#' \deqn{
#'   K(t) = -\log(\lambda - t) + \log(\lambda), \quad t < \lambda.
#' }
#'
#' By default, \code{ExponentialCGF} is vectorized for i.i.d. replicates. If you supply a
#' \code{tvec} of length \emph{m}, it sums the CGF across \emph{m} i.i.d. exponentials,
#' each with parameter \eqn{\lambda}.
#'
#' @format An object of class \code{CGF} (R6), with typical methods like
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' # Evaluate K, K1, K2 at t=0 for rate=2
#' #  -> K(0)= -log(2-0)+log(2)=  -log(2)+log(2)=0
#' ExponentialCGF$K(0, 2)
#'
#' #  -> K1(0)= 1/(2-0)= 0.5  (the mean for Exp with rate=2 is 1/2)
#' ExponentialCGF$K1(0, 2)
#'
#' #  -> K2(0)= 1/(2-0)^2= 0.25
#' ExponentialCGF$K2(0, 2)
#'
#' @export
ExponentialCGF <- createCGF_fromVectorisedFunctions(
  K_vectorized_func    = function(tvec, lambda) {
    # K(t)= -log(lambda - t) + log(lambda)
    # Summation across tvec for i.i.d. replicates
    -log(lambda - tvec) + log(lambda)
  },
  K1_vectorized_func   = function(tvec, lambda) { 1 / (lambda - tvec)  },
  K2_vectorized_func   = function(tvec, lambda) { 1 / (lambda - tvec)^2 },
  K3_vectorized_func   = function(tvec, lambda) { 2 / (lambda - tvec)^3 },
  K4_vectorized_func   = function(tvec, lambda) { 6 / (lambda - tvec)^4 },
  ineq_constraint_func = function(tvec, lambda) {
    # Enforce tvec < lambda => tvec - lambda < 0
    tvec - lambda
  },
  analytic_tvec_hat_func = function(y, lambda) {
    # 1/(lambda - t)= y => t= lambda - 1/y
    lambda - 1 / y
  },
  param_adaptor = function(x) x[1],
  op_name = "ExponentialCGF"
)










#' Create a Parametric Exponential CGF Object
#'
#' @description
#' Constructs a CGF object for the Exponential distribution where the rate parameter
#' \eqn{\lambda(\theta)} is derived from a user-supplied function. 
#'
#' @param rate A function(\code{theta}) -> numeric, returning a positive rate parameter \eqn{\lambda}.
#'  - Optionally, you can pass an **adaptor** object.
#'  - In many cases, an **adaptor** object can be used instead (especially if \eqn{\lambda} is constant or index-based).
#' @param iidReps A positive integer specifying the number of i.i.d. replicates 
#'   the CGF handles in \code{tvec}. Defaults to 1. If \code{iidReps>1}, \code{tvec} must be length \eqn{iidReps}.
#' @param ... Additional arguments passed to the base CGF creation function, such as optional method overrides.
#'
#'
#' @return An object of class \code{CGF} 
#'
#' @examples
#' rate_func <- function(theta) theta[1]  # e.g. param -> 2
#' expo_cgf <- ExponentialModelCGF(rate=rate_func, iidReps=1)
#' # Evaluate K1 at t=0 with param=c(2)
#' expo_cgf$K1(0, c(2)) 
#'
#' @export
ExponentialModelCGF <- function(rate, iidReps=1, ...) {
  # Validate that `rate` is a single-argument function (or adaptor)
  rate_func <- validate_function_or_adaptor(rate)
  if (!is.numeric(iidReps) || length(iidReps)!=1 || iidReps<1) stop("'iidReps' must be a positive integer.")
  
  
  # Helper to check tvec length
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) {
      stop(sprintf("`tvec` must have length %d, got %d.", iidReps, length(tvec)))
    }
  }
  
  createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      -log(param[1] - tvec) + log(param[1])
    },
    K1_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      1 / (param[1] - tvec)
    },
    K2_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      1 / (param[1] - tvec)^2
    },
    K3_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      2 / (param[1] - tvec)^3
    },
    K4_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      6 / (param[1] - tvec)^4
    },
    ineq_constraint_func = function(tvec, param) {
      # ensure tvec < lam => tvec - lam < 0
      tvec - param[1]
    },
    param_adaptor = rate_func,
    analytic_tvec_hat_func = function(y, param) {
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      param[1] - 1 / y
    },
    op_name = "ExponentialModelCGF",
    ...
  )
}