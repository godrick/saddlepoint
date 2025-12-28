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





#' @noRd
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
    # K2_elem     = function(tvec, pm) {
    #   num <- (1 - pm[,2]) * exp(tvec)
    #   denom <- 1 - (1 - pm[,2]) * exp(tvec)
    #   pm[,1] * num / denom^2
    # },
    K2_elem     = function(tvec, pm) {
      num <- (1 - pm[,2]) * exp(tvec)
      denom <- 1 - (1 - pm[,2]) * exp(tvec)
      denom2 <- denom * denom
      pm[,1] * num / denom2
    },
    # K3_elem     = function(tvec, pm) {
    #   num <- (1 - pm[,2]) * exp(tvec) * (1 + (1 - pm[,2]) * exp(tvec))
    #   denom <- 1 - (1 - pm[,2]) * exp(tvec)
    #   pm[,1] * num / denom^3
    # },
    K3_elem     = function(tvec, pm) {
      num <- (1 - pm[,2]) * exp(tvec) * (1 + (1 - pm[,2]) * exp(tvec))
      denom <- 1 - (1 - pm[,2]) * exp(tvec)
      denom2 <- denom * denom
      denom3 <- denom2 * denom
      pm[,1] * num / denom3
    },
    # K4_elem     = function(tvec, pm) {
    #   e_t <- exp(tvec)
    #   e_2t <- exp(2 * tvec)
    #   alpha <- 1 - (1 - pm[,2]) * e_t
    #   bracket <- 1 + e_2t + 4 * e_t - 2 * pm[,2] * e_2t - 4 * pm[,2] * e_t + pm[,2]^2 * e_2t
    #   num <- (1 - pm[,2]) * e_t * bracket
    #   pm[,1] * num / (alpha^4)
    # },
    K4_elem     = function(tvec, pm) {
      e_t <- exp(tvec)
      e_2t <- exp(2 * tvec)
      alpha <- 1 - (1 - pm[,2]) * e_t
      alpha2 <- alpha * alpha
      alpha4 <- alpha2 * alpha2
      p2 <- pm[,2] * pm[,2]
      bracket <- 1 + e_2t + 4 * e_t - 2 * pm[,2] * e_2t - 4 * pm[,2] * e_t + p2 * e_2t
      num <- (1 - pm[,2]) * e_t * bracket
      pm[,1] * num / alpha4
    },
    t_hat_elem  = function(x, pm) {
      q <- 1 - pm[,2]
      log(x) - log(q * (pm[,1] + x))
    },
    split_param_to_mat = split_rp,
    simulate_func = function(iidReps, parameter_vector, ...) {
      pm <- split_rp(parameter_vector)
      r <- as.numeric(pm[, 1])  # "size"
      p <- as.numeric(pm[, 2])  # "prob"

      if (any(!is.finite(r)) || any(r <= 0)) stop("NegBinCGF$rsim: 'r' (size) must be finite and > 0.")

      if (any(!is.finite(p)) || any(p <= 0 | p > 1)) stop("NegBinCGF$rsim: 'p' must be finite and in (0, 1].")


      d <- length(r)
      out <- stats::rnbinom(
        n    = d * iidReps,
        size = rep.int(r, times = iidReps),
        prob = rep.int(p, times = iidReps)
      )
      matrix(out, nrow = d, ncol = iidReps)
    },
    iidReps = iidReps,
    op_name = op_name,
    ineq_elem = function(tvec, pm) (1 - pm[,2]) * exp(tvec) - 1,
    ...
  )
}


#' Negative Binomial CGF Object
#'
#' A ready-to-use CGF object for the Negative Binomial distribution with
#' number of successes \eqn{r} and success probability \eqn{p}.
#' The \code{parameter_vector} is \eqn{c(r, p)}, and the actual
#' CGF is \deqn{K(t) = r\,[\log p - \log(1 - (1-p)\,e^t)], \quad t < -\log(1-p).}
#'
#'
#' @format An object of class \code{CGF} (R6) with methods \code{K}, \code{K1},
#' \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' NegBinCGF$K1(0, c(10, 0.25))  # E[X] = r(1-p)/p = 10 * 0.75 / 0.25 = 30
#'
#' @export
NegBinCGF <- .negbin_base_cgf(iidReps = "any", op_name = "NegBinCGF")



















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
  base_cgf <- .negbin_base_cgf(iidReps = iidReps, op_name = "NegBinModelCGF", ...)
  adaptCGF(
    cgf = base_cgf,
    adaptor = function(theta) { c(r_fn(theta), p_fn(theta)) }
  )
}















