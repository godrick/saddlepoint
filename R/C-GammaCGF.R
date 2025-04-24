# R/GammaCGF.R
# Objects: GammaCGF, GammaModelCGF

#' Gamma CGF Object
#'
#' A ready-to-use CGF object for the Gamma distribution with shape \eqn{\alpha}
#' and rate \eqn{\beta}. The \code{parameter_vector} used when calling methods such as `K(tvec, parameter_vector)`
#' should be a numeric vector \eqn{c(\alpha, \beta)}.
#' 
#'
#' @details
#' **CGF**: For a Gamma random variable \eqn{X} with shape \eqn{\alpha} and rate
#' \eqn{\beta}, the cumulant generating function is:
#' \deqn{K(t;\alpha, \beta) = -\alpha \,\log \bigl(1 - t/\beta\bigr), \quad t < \beta.}
#'
#' **Parameter Vector**: The \code{parameter_vector} is assumed to have the form 
#' \eqn{(\alpha, \beta)}. You must ensure 
#' that \eqn{t < \beta} for valid evaluations.
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
  ineq_constraint_func = function(tvec, param) {
    beta <- param[2]
    tvec - beta
  },
  
  analytic_tvec_hat_func = function(x, param)  { param[2] - (param[1] / x) },
  op_name       = "GammaCGF"
)

# #' @noRd
# validateGammaLengths <- function(vec, param, iidReps) {
#   d <- length(param) / 2
#   if (!is.null(iidReps)) {
#     expected_len <- d * iidReps
#     if (length(vec) != expected_len) {
#       stop(sprintf("Length of tvec/x is %d; expected %d (parameter dimension d = %d, iidReps = %s).",
#                    length(vec), expected_len, d, iidReps))
#     }
#   } else if (length(vec) %% d != 0) {
#     stop(sprintf("Length of tvec/x (%d) is not a multiple of the parameter dimension (%d).",
#                  length(vec), d))
#   }
# }



#' @noRd
.GammaModelCGF_internal <- function(iidReps, ...) {

  createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      -alpha * log1p(-tvec / beta)
    },
    K1_vectorized_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      alpha / (beta - tvec)
    },
    K2_vectorized_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      alpha / (beta - tvec)^2
    },
    K3_vectorized_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      2*alpha / (beta - tvec)^3
    },
    K4_vectorized_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      6*alpha / (beta - tvec)^4
    },
    ineq_constraint_func = function(tvec, sr) {
      validate2ParsLengths(tvec, sr, iidReps)
      len_s <- length(sr)/2
      beta  <- sr[len_s + 1:len_s]
      tvec - beta
    },
    analytic_tvec_hat_func = function(x, sr) {
      validate2ParsLengths(x, sr, iidReps)
      len_s <- length(sr)/2
      alpha <- sr[1:len_s]
      beta  <- sr[len_s + 1:len_s]
      beta - (alpha / x)
    },
    op_name = "GammaModelCGF",
    ...
  )
}





#' Create a Parametric Gamma CGF Object
#'
#' @description
#' Creates a CGF object for the Gamma distribution with shape \eqn{\alpha(\theta)} and 
#' rate \eqn{\beta(\theta)} defined by user-provided parameter functions. 
#' This function supports both i.i.d. and non-identical usage.
#' 
#'  
#'
#' @param shape A function (or `adaptor`)  that accepts a single parameter vector \code{theta} and returns the shape parameter.
#' @param rate A function (or `adaptor`)  that accepts a single parameter vector \code{theta} and returns the rate parameter.
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional arguments passed to the underlying CGF creation function.
#'
#'
#' @return A `CGF` object.
#' @export
GammaModelCGF <- function(shape,
                          rate,
                          iidReps = "any",
                          ...) {
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  
  shape_fn <- validate_function_or_adaptor(shape)
  rate_fn  <- validate_function_or_adaptor(rate)
  
  base_cgf <- .GammaModelCGF_internal(iidReps, ...)
  adaptCGF(
    cgf = base_cgf, 
    adaptor = function(theta) c(shape_fn(theta), rate_fn(theta))
  )
}



