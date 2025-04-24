# R/PoissonCGF.R
# Objects: PoissonCGF, validateLamLengths, .PoissonModelCGF_internal, PoissonModelCGF





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
  K4_vectorized_func = function(tvec, lambda) { lambda[1] * exp(tvec) }
  , analytic_tvec_hat_func = function(x, lambda)  log(x / lambda[1]) # solution of the saddlepoint equation
  , op_name = "PoissonCGF"
)





















# ------------------------------------------------------------------
# Helper: validate lam to match tvec length, 
# depending on the presence/value of iidReps
# ------------------------------------------------------------------

#' Validate Compatibility of Lambda and Vector Lengths
#'
#' @description
#' Validates that a numeric vector `lam` (Poisson rates) is 
#' compatible with the length of another numeric vector (`vec`, such as `tvec` or `x`)
#' based on the value of `iidReps`. It ensures:
#' - For an i.i.d. scenario (`iidReps` specified): 
#'   `length(vec)` equals `length(lam) * iidReps`.
#' - For a non-i.i.d. scenario (`iidReps` is NULL): 
#'   `length(vec)` is a multiple of `length(lam)`.
#'
#' @param vec Numeric vector (e.g., `tvec` or `x`) whose length is checked.
#' @param lam Numeric vector of Poisson rates.
#' @param iidReps Either `NULL` or a positive integer specifying replication.
#'
#' @return No return value; stops execution if validation fails.
#'
#' @noRd
validateLamLengths <- function(vec, lam, iidReps) {
  len_vec <- length(vec)
  len_lam <- length(lam)
  
  if (!is.null(iidReps)) {
    expected_ <- len_lam * iidReps
    if (len_vec != expected_) {
      stop(sprintf(
        "`(tvec/x)` has length %d; expected %d (lam length=%d, iidReps=%d).",
        len_vec, expected_, len_lam, iidReps
      ))
    }
  } else if (len_vec %% len_lam != 0) {
    stop(sprintf(
      "Length mismatch: lambda length=%d, (tvec/x) length=%d. Length of (tvec/x) must be a multiple of lambda length.",
      len_lam, len_vec
    ))
  }
  # No return value needed; function ends if validation passes.
}












.PoissonModelCGF_internal <- function(iidReps, ...) {
  
  # Create the base Poisson CGF using pre-validated iidReps
  createCGF_fromVectorisedFunctions(
    K_vectorized_func  = function(tvec, lam) { 
      validateLamLengths(tvec, lam, iidReps)
      lam*(exp(tvec) - 1)
    },
    K1_vectorized_func = function(tvec, lam) { 
      validateLamLengths(tvec, lam, iidReps)
      lam*exp(tvec) 
    },
    K2_vectorized_func = function(tvec, lam) { 
      validateLamLengths(tvec, lam, iidReps)
      lam*exp(tvec) 
    },
    K3_vectorized_func = function(tvec, lam) { 
      validateLamLengths(tvec, lam, iidReps)
      lam*exp(tvec) 
    },
    K4_vectorized_func = function(tvec, lam) { 
      validateLamLengths(tvec, lam, iidReps)
      lam*exp(tvec) 
    },
    analytic_tvec_hat_func = function(x, lam)  {
      validateLamLengths(x, lam, iidReps)
      log(x / lam)
    }, 
    op_name = "PoissonModelCGF",
    ...  # Pass any additional optional arguments
  )
  
}













#' Create a Parametric Poisson CGF Object
#'
#' @description
#' Creates a CGF for the Poisson distribution with a rate parameter (\eqn{\lambda}) determined 
#' by a user-supplied function or adaptor. Supports i.i.d. and non-identical contexts 
#' with optional length enforcement via `iidReps`.
#' 
#' @param lambda An `adaptor` or a function mapping a parameter vector to a numeric vector rate(s).
#'    The function must accept a single parameter vector and return the Poisson rate parameter (\eqn{\lambda}) or vector of rates.
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
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
#' # ex 2: non-identical scenario with replication: lambda returns c(2,5), iidReps=3
#' ex_repeat <- PoissonModelCGF(lambda = adaptor(indices = 1:2), iidReps = 3)
#' # OR ex_repeat <- PoissonModelCGF(lambda = function(x) c(2,5) )
#' ex_repeat$K1(rep(0,6), c(2,5))
#' 
#' @export
PoissonModelCGF <- function(lambda, iidReps = "any", ...) {
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  
  # Validate lambda function/adaptor
  lambda_adaptor <- validate_function_or_adaptor(obj = lambda)
  
  base_cgf <- .PoissonModelCGF_internal(iidReps, ...)
  
  adaptCGF(cgf = base_cgf, adaptor = lambda_adaptor)
}


