# R/C-NegBinCGF.R



# # G(z) = ( p / [1 - (1-p) z] )^r
# NB_pgf <- function(z, r, p) {
#   (p / (1 - (1-p)*z))^r
# }
# 
# # d/dz [ G(z) ] = G(z) * [ r (1-p) / (1 - (1-p)*z ) ]
# NB_pgf_prime <- function(z, r, p) {
#   r*(1-p)*NB_pgf(z, r, p) / (1 - (1-p)*z)
# }




#' Negative Binomial CGF Object
#'
#' A ready-to-use CGF object for the Negative Binomial distribution with number of successes \eqn{r}
#' and success probability parameter \eqn{p}. The \code{parameter_vector} used when calling methods such as `K(tvec, parameter_vector)`
#' should be a numeric vector \eqn{c(r, p)}.
#' 
#'
#' @format An object of class \code{CGF} (an R6 class), with the usual methods:
#' \code{K, K1, K2, K3operator, K4operator}, etc.
#'
#'
#' @export
NegBinCGF <- createCGF_fromVectorisedFunctions(
  
  # ------------------------------------------------------------------------
  # K(tvec, r, p):
  # K(t) = r * log( p / [1 - (1-p)*exp(t)] )
  # ------------------------------------------------------------------------
  K_vectorized_func = function(tvec, param) {
    r <- param[1]
    p <- param[2]
    # alpha(t) = 1 - (1-p)*exp(t)
    # => K(t) = r * [ log(p) - log(alpha(t)) ]
    alpha <- 1 - (1 - p)*exp(tvec)
    r*(log(p) - log(alpha))
  },
  K1_vectorized_func = function(tvec, param) {
    r <- param[1]
    p <- param[2]
    # K'(tvec, r, p):
    #   = r * (1-p)*exp(t) / [alpha(t)]
    num_ <- (1-p)*exp(tvec)
    denom_ <- 1 - (1 - p)*exp(tvec)
    r*num_/denom_
  },
  K2_vectorized_func = function(tvec, param) {
    r <- param[1]
    p <- param[2]
    # K''(tvec, r, p):
    #   = r*(1-p)*exp(t) / [alpha(t)]^2
    num_ <- (1-p)*exp(tvec)
    denom_ <- 1 - (1 - p)*exp(tvec)
    r*num_/denom_^2
  },
  K3_vectorized_func = function(tvec, param) {
    r <- param[1]
    p <- param[2]
    # K'''(tvec, r, p):
    # = r*(1-p)*exp(t) * [1 + (1-p)*exp(t)] / [alpha(t)]^3
    num_ <- (1-p)*exp(tvec) * (1 + (1-p)*exp(tvec))
    denom_ <- 1 - (1 - p)*exp(tvec)
    r*num_/denom_^3
  },
  K4_vectorized_func = function(tvec, param) {
    r <- param[1]
    p <- param[2]
    # K''''(tvec, r, p):
    #   = r*(1-p)*exp(t) * [1 + exp(2t) + 4exp(t)
    #                       - 2p exp(2t) - 4p exp(t) + p^2 exp(2t)] / [alpha]^4
    e_t <- exp(tvec)
    e_2t <- exp(2*tvec)
    alpha_ <- 1 - (1 - p)*e_t
    
    bracket_ <- 1 + e_2t + 4*e_t - 2*p*e_2t - 4*p*e_t + p^2*e_2t
    numerator_ <- (1 - p)*e_t * bracket_
    r*numerator_/(alpha_^4)
  },
  
  # We want tvec < beta => tvec - beta < 0
  # We'll just return (tvec - beta), which must be negative
  ineq_constraint_func = function(tvec, param) {
    # The domain is:  (1-p)*exp(t) < 1.
    (1 - param[2])*exp(tvec) - 1
  },
  
  analytic_tvec_hat_func = function(x, param)  { 
    # solve K'(t_i)= x_i => r q_i e^t / [1 - q_i e^t ]= x_i
    # => e^t= x_i / [q_i (r_i + x_i)]
    q <- 1 - param[2]
    log(x) - log(q*(param[1] + x))
  },
  op_name       = "NegBinCGF"
)







#' @noRd
validate2ParsLengths <- function(vec, param, iidReps) {
  d <- length(param) / 2
  if (!is.null(iidReps)) {
    expected_len <- d * iidReps
    if (length(vec) != expected_len) {
      stop(sprintf("Length of tvec/x is %d; expected %d (parameter dimension d = %d, iidReps = %s).",
                   length(vec), expected_len, d, iidReps))
    }
  } else if (length(vec) %% d != 0) {
    stop(sprintf("Length of tvec/x (%d) is not a multiple of the parameter dimension (%d).",
                 length(vec), d))
  }
}






#' @noRd
.NegBinModelCGF_internal <- function(iidReps, ...) {
  createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      
      alpha <- 1 - (1 - p)*exp(tvec)
      r*(log(p) - log(alpha))
    },
    K1_vectorized_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      
      num_ <- (1-p)*exp(tvec)
      denom_ <- 1 - (1 - p)*exp(tvec)
      r*num_/denom_
    },
    K2_vectorized_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      
      num_ <- (1-p)*exp(tvec)
      denom_ <- 1 - (1 - p)*exp(tvec)
      r*num_/denom_^2
    },
    K3_vectorized_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      
      num_ <- (1-p)*exp(tvec) * (1 + (1-p)*exp(tvec))
      denom_ <- 1 - (1 - p)*exp(tvec)
      r*num_/denom_^3
    },
    K4_vectorized_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      
      e_t <- exp(tvec)
      e_2t <- exp(2*tvec)
      alpha_ <- 1 - (1 - p)*e_t
      
      bracket_ <- 1 + e_2t + 4*e_t - 2*p*e_2t - 4*p*e_t + p^2*e_2t
      numerator_ <- (1 - p)*e_t * bracket_
      r*numerator_/(alpha_^4)
    },
    

    ineq_constraint_func = function(tvec, param) {
      validate2ParsLengths(tvec, param, iidReps)
      len_r <- length(param)/2
      p <- param[len_r + 1:len_r]
      (1 - p)*exp(tvec) - 1
    },
    
    analytic_tvec_hat_func = function(x, param)  { 
      validate2ParsLengths(x, param, iidReps)
      len_r <- length(param)/2
      r <- param[1:len_r]
      p <- param[len_r + 1:len_r]
      q <- 1 - p
      log(x) - log(q*(r + x))
    },
    op_name       = "NegBinModelCGF"
  )
}








#' Create a Parametric Negative Binomial CGF Object
#'
#' @description
#' Creates a CGF object for the Negative Binomial distribution with number of successes \eqn{r(\theta)} and 
#' success probability parameter \eqn{p(\theta)} defined by user-provided parameter functions. 
#' This function supports both i.i.d. and non-identical usage.
#' 
#'  
#'
#' @param r A function (or `adaptor`)  that accepts a single parameter vector \code{theta} and returns the number of successes.
#' @param p A function (or `adaptor`)  that accepts a single parameter vector \code{theta} and returns a scalar success probability or a vector of success probabilities.
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional arguments passed to the underlying CGF creation function.
#'
#'
#' @return A `CGF` object.
#' @export
NegBinModelCGF <- function(r, p, iidReps = "any", ...) {
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  r_fn <- validate_function_or_adaptor(r)
  p_fn <- validate_function_or_adaptor(p)
  base_cgf <- .NegBinModelCGF_internal(iidReps, ...)
  adaptCGF(
    cgf = base_cgf,
    adaptor = function(theta) { c(r_fn(theta), p_fn(theta)) }
  )
}
















