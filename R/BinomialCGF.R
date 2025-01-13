# R/BinomialCGF.R

#' Binomial CGF Object
#'
#' A ready-to-use CGF object for the Binomial distribution.
#' 
#' @details
#' **Parameter Vector**:
#' When calling methods like \code{K}, \code{K1}, etc. on \code{BinomialCGF}, you must provide
#' a numeric vector \code{param = c(n, p)}, where:
#' \itemize{
#'   \item \eqn{n} is the number of trials (must be a positive integer),
#'   \item \eqn{p} is the probability of success (must be in (0,1)).
#' }
#' 
#' By default, \code{BinomialCGF} supports vectorized evaluation for i.i.d. replicates of
#' Binomial random variables. If you provide a \code{tvec} of length \eqn{m}, it will
#' compute the sum of CGFs for \eqn{m} independent copies, each with parameters \eqn{(n, p)}.
#' 
#' 
#' @format An object of class \code{CGF} (R6), with standard methods 
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.
#'
#' @examples
#' # Evaluate K1, K2 at t = 0 for n = 10, p = 0.5
#' # This should return the expected value and variance of a Binomial(10, 0.5) random variable.
#' BinomialCGF$K1(0, c(10, 0.5))
#' BinomialCGF$K2(0, c(10, 0.5))
#'
#' @export
BinomialCGF <- createCGF_fromVectorisedFunctions(
  K_vectorized_func = function(tvec, params) { 
    n <- params[1]
    p <- params[2]
    n*log1p(p*expm1(tvec))
  },
  K1_vectorized_func = function(tvec, params) { 
    n <- params[1]
    p <- params[2]
    n*p*exp(tvec) / (1 + p*expm1(tvec))
  },
  K2_vectorized_func = function(tvec, params) { 
    n <- params[1]
    p <- params[2]
    n*p*(1 - p)*exp(tvec) / (1 - p + p*exp(tvec))^2
  },
  K3_vectorized_func = function(tvec, params) { 
    n <- params[1]
    p <- params[2]
    tmp0 <- 1 - p + p*exp(tvec)
    n*(1 - p)*p*exp(tvec)*(1 - p - p*exp(tvec)) / (tmp0^3)
  },
  K4_vectorized_func = function(tvec, params) { 
    n <- params[1]
    p <- params[2]
    tmp0 <- p*exp(tvec)
    tmp1 <- 1-p
    tmp2 <- tmp1 + p*exp(tvec)
    n*tmp1*tmp0*(4*p*tmp0 + tmp0^2 - 4*tmp0 + tmp1^2) / (tmp2^4)
  },
  analytic_tvec_hat_func = function(y, params) { 
    # Solve K'(t) = y => t = log((y*(1-p)) / (p*(n-y)))
    n <- params[1]
    p <- params[2]
    # numeric check omitted for brevity, but typically you'd ensure y<n
    log(y*(1 - p)) - log(p*(n - y))
  },
  param_adaptor = function(x) x[1:2],
  op_name = "BinomialCGF"
)




#' Create a Parametric Binomial CGF Object
#'
#' #' @description
#' Constructs a Binomial CGF object where the number of trials \eqn{n(\theta)} and
#' the probability of success \eqn{p(\theta)} are derived from user-supplied logic. 
#' Each of \code{n} and \code{p} can be an R function (or an adaptor) that 
#' accepts a single parameter vector \code{theta} and returns the respective value.
#' 
#' @param n A function(\code{theta}) returning a **positive integer** \eqn{n} (the number of trials).
#'   - Typically, you define a simple function, e.g. \code{function(theta) theta[1]}.
#'   - In many cases, an **adaptor** object can be used instead (especially if \eqn{n} is constant or index-based).
#' @param p A function(\code{theta}) returning a numeric \eqn{p} in \eqn{[0, 1]}.
#'   - Similarly, can be a direct R function or an adaptor object for constant or index-based usage.
#' @param iidReps A positive integer. Specifies the number of i.i.d. replicates the CGF will expect in \code{tvec}.
#' @param ... Additional arguments passed to the base CGF creation function, such as optional method overrides.
#'
#'
#' @return A CGF object
#' 
#' @examples
#' # Simple usage:
#' n_func <- function(theta) theta[1]  # e.g., 10
#' p_func <- function(theta) theta[2]  # e.g., 0.3
#' cgf_obj <- BinomialModelCGF(n = n_func, p = p_func, iidReps=1)
#' # Evaluate K1 at t=0 with param=c(10,0.3)
#' cgf_obj$K1(0.5, c(10, 0.3))
#'
#' # Using adaptor objects 
#' n_adapt <- adaptor(fixed_param=10)  # always 10
#' p_adapt <- adaptor(indices=1)       # interpret param[1] as p
#' cgf_obj2 <- BinomialModelCGF(n=n_adapt, p=p_adapt, iidReps=1)
#' cgf_obj2$K1(0, 0.3) # p=0.3 from param
#' 
#' @export
BinomialModelCGF <- function(n, p, iidReps = 1, ...) {
  n_fn <- validate_function_or_adaptor(n)
  p_fn <- validate_function_or_adaptor(p)
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1 || iidReps != as.integer(iidReps)) stop("'iidReps' must be a positive integer.")
  
  # Helper function to validate tvec length
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) stop(sprintf("`tvec` must have length %d (got %d).", iidReps, length(tvec)))
  }
  
  createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, params) {
      check_tvec(tvec)
      n_val <- params[1]
      p_val <- params[2]
      n_val*log1p(p_val*expm1(tvec))
    },
    K1_vectorized_func = function(tvec, params) {
      check_tvec(tvec)
      n_val <- params[1]
      p_val <- params[2]
      n_val * p_val * exp(tvec) / (1 + p_val * expm1(tvec))
    },
    K2_vectorized_func = function(tvec, params) {
      check_tvec(tvec)
      n_val <- params[1]
      p_val <- params[2]
      n_val * p_val * (1 - p_val) * exp(tvec) / (1 + p_val * (exp(tvec) - 1))^2
    },
    K3_vectorized_func = function(tvec, params) {
      check_tvec(tvec)
      n_val <- params[1]
      p_val <- params[2]
      tmp0  <- 1 - p_val + p_val*exp(tvec)
      n_val*(1 - p_val)*p_val*exp(tvec)*(1 - p_val - p_val*exp(tvec)) / (tmp0^3)
    },
    K4_vectorized_func = function(tvec, params) {
      check_tvec(tvec)
      n_val <- params[1]
      p_val <- params[2]
      tmp0  <- p_val*exp(tvec)
      tmp1  <- 1-p_val
      tmp2  <- tmp1 + p_val*exp(tvec)
      n_val*tmp1*tmp0*(4*p_val*tmp0 + tmp0^2 - 4*tmp0 + tmp1^2) / (tmp2^4)
    },
    param_adaptor = function(theta) c(n_fn(theta), p_fn(theta)),
    analytic_tvec_hat_func = function(y, params) { 
      n_val <- params[1]
      p_val <- params[2]
      # Compute t = log((y * (1 - p)) / (p * (n - y)))
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      log(y * (1 - p_val)) - log(p_val * (n_val - y))
    },
    op_name = "BinomialModelCGF",
    ...
  )
}



