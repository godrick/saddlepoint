# R/PoissonCGF.R
# Objects: PoissonCGF, PoissonModelCGF, expandLam





#' Poisson CGF Object
#'
#' A ready-to-use CGF object for the Poisson distribution. 
#' 
#' #' @details
#' **Parameter Vector**: The `parameter_vector` for \code{PoissonCGF} is 
#' interpreted as the rate \eqn{\lambda}. By default, this object is vectorized 
#' for i.i.d. replicates of Poisson variables.
#' 
#' @format An object of class \code{CGF} (an R6 class), with standard methods 
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.
#'
#' @examples
#' # Evaluate K at t = 0.1 for lambda = 2
#' # PoissonCGF$K(0.1, 2)
#' 
#' @export
PoissonCGF <- createCGF_fromVectorisedFunctions(
  K_vectorized_func = function(tvec, lambda) { lambda[1] * (exp(tvec) - 1) },
  K1_vectorized_func = function(tvec, lambda) { lambda[1] * exp(tvec) },
  K2_vectorized_func = function(tvec, lambda) { lambda[1] * exp(tvec) },
  K3_vectorized_func = function(tvec, lambda) { lambda[1] * exp(tvec) },
  K4_vectorized_func = function(tvec, lambda) { lambda[1] * exp(tvec) },
  ,analytic_tvec_hat_func = function(y, lambda)  log(y / lambda[1]) # solution of the saddlepoint equation
  , op_name = "PoissonCGF"
)














# ------------------------------------------------------------------
# Helper: expand lam to match tvec length, 
# depending on the presence/value of iidReps
# ------------------------------------------------------------------

#' Expand a Poisson rate vector to match vector length
#'
#' @description
#' Expands a numeric vector `lam` (Poisson rates) to match the length of another
#' numeric vector (`tvec` or `y`), based on the value of `iidReps`. It handles:
#' - i.i.d. scenario when `iidReps` is specified.
#' - Non-identical/flexible scenario when `iidReps` is NULL.
#'
#' @param vec Numeric vector (e.g., `tvec` or `y`) whose length determines the expansion.
#' @param lam Numeric vector of Poisson rates.
#' @param iidReps Either `NULL` or a positive integer specifying replication.
#'
#' @return A numeric vector of Poisson rates expanded to match the length of `vec`.
#'
#' @noRd
expandLam <- function(vec, lam, iidReps) {
  len_vec <- length(vec)
  len_lam <- length(lam)
  
  if (!is.null(iidReps)) {
    expected_ <- len_lam * iidReps
    if (len_vec != expected_) {
      stop(sprintf(
        "`tvec/observations` has length %d; expected %d (lam length=%d, iidReps=%d).",
        len_vec, expected_, len_lam, iidReps
      ))
    }
    return( rep(lam, times = iidReps) ) # This covers both univariate and multivariate i.i.d. cases when we want to replicate the lam vector as blocks.
  }
  
  if (is.null(iidReps) && len_vec %% len_lam == 0){
    times_ <- len_vec / len_lam
    return( rep(lam, times = times_) )
  }
      
  stop(sprintf(
    "Length mismatch: lambda length=%d, tvec/observation length=%d. Either lambda=1, tvec/observation length is a multiple of lambda length, or they match exactly.",
    len_lam, len_vec
  ))
  
}



























#' Create a Parametric Poisson CGF Object
#'
#' @description
#' Creates a CGF for the Poisson distribution with a rate parameter (\eqn{\lambda}) determined 
#' by a user-supplied `lambda` function or adaptor. Supports i.i.d. and non-identical contexts 
#' with optional length enforcement via `iidReps`.
#' 
#' @param lambda An `adaptor` or a function mapping a parameter vector to a numeric vector rate(s).
#'    The function must accept a single parameter vector and return the Poisson rate parameter (\eqn{\lambda}) or vector of rates.
#' @param iidReps Either \code{NULL} (no length restriction on \code{tvec}) or a positive integer
#'   specifying the required length of \code{tvec}. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to the base CGF creation function, such as optional method overrides.
#'
#' @return A CGF object.
#' 
#' 
#' @examples
#' # ex 1: i.i.d. scenario for univariate Poisson(rate = 2)
#' ex_iid <- PoissonModelCGF(lambda = function(x) x[1])
#' # OR ex_iid <- PoissonModelCGF(lambda = function(x) x[1], iidReps = 3)
#' ex_iid$K1(rep(0,3), 2)
#'
#' # ex 3: non-identical scenario with replication: lambda returns c(2,5), iidReps=3
#' ex_repeat <- PoissonModelCGF(lambda = adaptor(indices = 1:2), iidReps = 3)
#' # OR ex_repeat <- PoissonModelCGF(lambda = function(x) c(2,5) )
#' ex_repeat$K1(rep(0,6), c(2,5))
#' 
#' @export
PoissonModelCGF <- function(lambda, iidReps = NULL, ...) {
  # validate iidReps if not NULL
  if (!is.null(iidReps) && (!is.numeric(iidReps) || length(iidReps) != 1 ||
                            iidReps < 1 || iidReps != as.integer(iidReps))) {
    stop("`iidReps` must be NULL or a positive integer.")
  }
  
  # Validate lambda function/adaptor
  lambda_adaptor <- validate_function_or_adaptor(obj = lambda)
  
  
  # Create the base Poisson CGF
  base_cgf <- createCGF_fromVectorisedFunctions(
    K_vectorized_func  = function(tvec, lam) { 
      lam_expanded <- expandLam(tvec, lam, iidReps)
      lam_expanded*(exp(tvec) - 1)
    },
    K1_vectorized_func = function(tvec, lam) { 
      lam_expanded <- expandLam(tvec, lam, iidReps)
      lam_expanded*exp(tvec) 
    },
    K2_vectorized_func = function(tvec, lam) { 
      lam_expanded <- expandLam(tvec, lam, iidReps)
      lam_expanded*exp(tvec) 
    },
    K3_vectorized_func = function(tvec, lam) { 
      lam_expanded <- expandLam(tvec, lam, iidReps)
      lam_expanded*exp(tvec) 
    },
    K4_vectorized_func = function(tvec, lam) { 
      lam_expanded <- expandLam(tvec, lam, iidReps)
      lam_expanded*exp(tvec) 
    },
    analytic_tvec_hat_func = function(y, lam)  {
      lam_expanded <- expandLam(y, lam, iidReps)
      log(y / lam_expanded)
    }, 
    op_name = "PoissonModelCGF",
    ...  # Pass any additional optional arguments
  )
  adaptCGF(cgf = base_cgf, param_adaptor = lambda_adaptor)
}


