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
## Internal factory – build NegBin CGF via univariate utilities
.negbin_base_cgf <- function(iidReps, op_name, ...) {
  split_rp <- function(param) {
    ln <- length(param) / 2
    if (ln != as.integer(ln) || ln < 1) stop("NegBin params must be concatenated as c(r[1:L], p[1:L]).")
    idx <- seq_len(ln)
    cbind(r = param[idx], p = param[ln + idx])
  }

  .make_univariate_model_cgf_matrix(
    K_elem      = function(tvec, pm) {
      alpha <- 1 - (1 - pm[,2]) * exp(tvec)
      pm[,1] * (log(pm[,2]) - log(alpha))
    },
    K1_elem     = function(tvec, pm) {
      num <- (1 - pm[,2]) * exp(tvec)
      denom <- 1 - (1 - pm[,2]) * exp(tvec)
      pm[,1] * num / denom
    },
    K2_elem     = function(tvec, pm) {
      num <- (1 - pm[,2]) * exp(tvec)
      denom <- 1 - (1 - pm[,2]) * exp(tvec)
      pm[,1] * num / denom^2
    },
    K3_elem     = function(tvec, pm) {
      num <- (1 - pm[,2]) * exp(tvec) * (1 + (1 - pm[,2]) * exp(tvec))
      denom <- 1 - (1 - pm[,2]) * exp(tvec)
      pm[,1] * num / denom^3
    },
    K4_elem     = function(tvec, pm) {
      e_t <- exp(tvec)
      e_2t <- exp(2 * tvec)
      alpha <- 1 - (1 - pm[,2]) * e_t
      bracket <- 1 + e_2t + 4 * e_t - 2 * pm[,2] * e_2t - 4 * pm[,2] * e_t + pm[,2]^2 * e_2t
      num <- (1 - pm[,2]) * e_t * bracket
      pm[,1] * num / (alpha^4)
    },
    t_hat_elem  = function(x, pm) {
      q <- 1 - pm[,2]
      log(x) - log(q * (pm[,1] + x))
    },
    split_param_to_mat = split_rp,
    iidReps = iidReps,
    op_name = op_name,
    ineq_elem = function(tvec, pm) (1 - pm[,2]) * exp(tvec) - 1,
    ...
  )
}

#' @export
NegBinCGF <- .negbin_base_cgf(iidReps = "any", op_name = "NegBinCGF")







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
  .negbin_base_cgf(iidReps = iidReps, op_name = "NegBinModelCGF", ...)
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
  .check_iidReps(iidReps)
  r_fn <- validate_function_or_adaptor(r)
  p_fn <- validate_function_or_adaptor(p)
  base_cgf <- .NegBinModelCGF_internal(iidReps, ...)
  adaptCGF(
    cgf = base_cgf,
    adaptor = function(theta) { c(r_fn(theta), p_fn(theta)) }
  )
}















