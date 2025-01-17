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


#' Geometric CGF Object
#'
#' A ready-to-use CGF object for a single-parameter Geometric distribution.
#' This object correspond to the random variable that counts the number of failures before achieving the first success.
#' The `parameter_vector` for `GeometricCGF` is taken to be the success probability \eqn{prob}.
#' By default, this object is vectorized for i.i.d. replicates of Geometric random variables.
#'
#' @seealso \code{\link{GeometricModelCGF}}
#'
#' @format An object of class \code{CGF} (R6), with usual methods:
#' \code{K, K1, K2, K3operator, K4operator}, etc.
#'
#' @examples
#' # the expected value of X~Geometric(prob = 0.3) using the CGF 
#' GeometricCGF$K1(0, 0.3)
#'
#' @export
GeometricCGF <- createCGF_fromVectorisedFunctions(
  
  K_vectorized_func = function(tvec, p) log(p[1])-log( 1 - exp(tvec) + p[1]*exp(tvec) ),
  K1_vectorized_func = function(tvec, p) (exp(tvec) - p[1]*exp(tvec)) / (1 - exp(tvec) + p[1]*exp(tvec)),
  
  K2_vectorized_func = function(tvec, p) {
    tmp_ <- 1 - exp(tvec) + p[1]*exp(tvec)
    (exp(tvec) - p[1]*exp(tvec)) / tmp_^2 
  },
  
  K3_vectorized_func = function(tvec, p) {
    tmp_ <- 1 - exp(tvec) + p[1]*exp(tvec)
    (exp(tvec) - p[1]*exp(tvec)) * (1 + exp(tvec) - p[1]*exp(tvec) ) / tmp_^3
  },
  
  # K4(t; p) = as in snippet
  K4_vectorized_func = function(tvec, p) {
    tmp_ <- 1 - exp(tvec) + p[1]*exp(tvec)
    (exp(tvec) - p[1]*exp(tvec)) * (1 + exp(2*tvec) + 4*exp(tvec) - 2*p[1]*exp(2*tvec) - 4*p[1]*exp(tvec) + p[1]^2*exp(2*tvec)) / tmp_^4
  },
  
  ineq_constraint_func = function(tvec, p) {
    # tvec < -log(1-p)
    # The value of this function is constrained to be non-positive
    # tvec + log(1 - p) < 0
    # return tvec + log(1 - p)
    (1-p[1])*exp(tvec) - 1
  },
  analytic_tvec_hat_func = function(y, p) log(y) - log(1 + y - p[1] - p[1]*y),
  op_name = "GeometricCGF"
)











# ----------------------------------------------------------------------------
# Next, A parametric "Model", for i.i.d. and non-identical usage.
# ----------------------------------------------------------------------------

# Helper to expand prob if needed
#' @noRd
expandP <- function(vec, pval, iidReps) {
  len_vec <- length(vec)
  len_p <- length(pval)
  
  if (!is.null(iidReps)) {
    # forced dimension => expected = len_p * iidReps
    expected_ <- len_p * iidReps
    if (len_vec != expected_) {
      stop(sprintf(
        "`tvec/observations` has length %d; expected %d (p length=%d, iidReps=%d).",
        len_vec, expected_, len_p, iidReps
      ))
    }
    return(rep(pval, times = iidReps))
  }
  
  if (is.null(iidReps) && (len_vec %% len_p == 0)) {
    times_ <- len_vec / len_p
    return( rep(pval, times = times_) )
  }
  
  stop(sprintf(
    "Length mismatch: prob length=%d, tvec/observations length=%d. Either lambda=1, tvec/observatiosn length is a multiple of prob length, or they match exactly.",
    len_p, len_vec
  ))
  
  
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
#' @param iidReps Either `NULL` (no forced dimension) or a positive integer specifying
#'   that \code{tvec} is partitioned into `iidReps` blocks. Each block is one copy
#'   of the geometric variables (or multiple if `prob` is vector).
#' @param ... Additional arguments passed to the underlying CGF creation function
#'   (e.g., optional operator overrides).
#'
#' @return A CGF object
#' @export
GeometricModelCGF <- function(prob, iidReps = NULL, ...) {
  
  # Validate
  if (!is.null(iidReps)) {
    if (!is.numeric(iidReps) || length(iidReps) != 1 ||
        iidReps < 1 || iidReps != as.integer(iidReps)) {
      stop("'iidReps' must be NULL or a positive integer.")
    }
  }
  
  p_adaptor <- validate_function_or_adaptor(prob)
  
  
  
  # Build the base CGF for vectorized usage
  base_cgf <- createCGF_fromVectorisedFunctions(
    
    # K(t)= log p - log(1 - e^t + p e^t), but we must expand prob if it's multi
    K_vectorized_func = function(tvec, pval) {
      p_expanded <- expandP(tvec, pval, iidReps)
      log(p_expanded) - log(1 - exp(tvec) + p_expanded * exp(tvec))
    },
    
    K1_vectorized_func = function(tvec, pval) {
      p_expanded <- expandP(tvec, pval, iidReps)
      tmp_ <- 1 - exp(tvec) + p_expanded*exp(tvec)
      (exp(tvec) - p_expanded*exp(tvec)) / tmp_^2 
    },
    
    K2_vectorized_func = function(tvec, pval) {
      p_expanded <- expandP(tvec, pval, iidReps)
      tmp_ <- 1 - exp(tvec) + p_expanded*exp(tvec)
      (exp(tvec) - p_expanded*exp(tvec)) / tmp_^2 
    },
    
    K3_vectorized_func = function(tvec, pval) {
      p_expanded <- expandP(tvec, pval, iidReps)
      tmp_ <- 1 - exp(tvec) + p_expanded*exp(tvec)
      (exp(tvec) - p_expanded*exp(tvec)) * (1 + exp(tvec) - p_expanded*exp(tvec) ) / tmp_^3
    },
    
    K4_vectorized_func = function(tvec, pval) {
      p_expanded <- expandP(tvec, pval)
      tmp_ <- 1 - exp(tvec) + p_expanded*exp(tvec)
      (exp(tvec) - p_expanded*exp(tvec)) * (1 + exp(2*tvec) + 4*exp(tvec) - 2*p_expanded*exp(2*tvec) - 4*p_expanded*exp(tvec) + p_expanded^2*exp(2*tvec)) / tmp_^4
    },
    

    ineq_constraint_func = function(tvec, pval) {
      # We'll return block of length(tvec), requiring each entry t_i < -log(1 - pval_i)
      # That is: t_i + log(1 - pval_i) < 0
      p_expanded <- expandP(tvec, pval)
      (1-p_expanded)*exp(tvec) - 1  # must be <0
    },
    
    analytic_tvec_hat_func = function(y, pval) {
      p_expanded <- expandP(y, pval)
      log(y) - log(1 + y - p_expanded - p_expanded*y)
    },
    
    op_name = "GeometricModelCGF",
    ...
  )
  adaptCGF(cgf = base_cgf, param_adaptor = p_adaptor)
}

