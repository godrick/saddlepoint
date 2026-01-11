# R/GammaCGF.R
# Objects: GammaCGF, GammaModelCGF


#' @noRd
.gamma_base_cgf <- function(iidReps, op_name, ...) {
  split_ab <- function(param) {
    ln <- length(param) / 2
    if (ln != as.integer(ln) || ln < 1) stop("Gamma params must be concatenated as c(shape[1:L], rate[1:L]).")
    idx <- seq_len(ln)
    cbind(shape = param[idx], rate = param[ln + idx])
  }

  .make_univariate_model_cgf_matrix(
    K_elem      = function(tvec, pm) -pm[,1] * log1p(-tvec / pm[,2]),
    K1_elem     = function(tvec, pm)  pm[,1] / (pm[,2] - tvec),
    # K2_elem     = function(tvec, pm)  pm[,1] / (pm[,2] - tvec)^2,
    K2_elem     = function(tvec, pm)  {
      denom <- pm[,2] - tvec
      denom2 <- denom * denom
      pm[,1] / denom2
    },
    # K3_elem     = function(tvec, pm)  2 * pm[,1] / (pm[,2] - tvec)^3,
    K3_elem     = function(tvec, pm)  {
      denom <- pm[,2] - tvec
      denom2 <- denom * denom
      denom3 <- denom2 * denom
      2 * pm[,1] / denom3
    },
    # K4_elem     = function(tvec, pm)  6 * pm[,1] / (pm[,2] - tvec)^4,
    K4_elem     = function(tvec, pm)  {
      denom <- pm[,2] - tvec
      denom2 <- denom * denom
      denom4 <- denom2 * denom2
      6 * pm[,1] / denom4
    },
    t_hat_elem  = function(x, pm) pm[,2] - (pm[,1] / x),
    split_param_to_mat = split_ab,
    rsim_elem = function(n, tvec, pm, ...) {
      shape <- as.numeric(pm[, 1])
      rate  <- as.numeric(pm[, 2])

      if (any(!is.finite(shape)) || any(shape <= 0)) {
        stop("GammaCGF$rsim: 'shape' must be finite and > 0.")
      }
      if (any(!is.finite(rate)) || any(rate <= 0)) {
        stop("GammaCGF$rsim: 'rate' must be finite and > 0.")
      }

      rate_tilt <- rate - tvec
      if (any(!is.finite(rate_tilt)) || any(rate_tilt <= 0)) {
        stop("GammaCGF$rsim: 'rate - tvec' must be finite and > 0.", call. = FALSE)
      }

      vector_length <- length(tvec)
      matrix(
        stats::rgamma(
          n = n * vector_length,
          shape = rep.int(shape, times = n),
          rate = rep.int(rate_tilt, times = n)
        ),
        nrow = vector_length,
        ncol = n
      )
    },

    iidReps = iidReps,
    op_name = op_name,
    ineq_elem = function(tvec, pm) tvec - pm[,2],
    ...
  )
}


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
GammaCGF <- .gamma_base_cgf(iidReps = "any", op_name = "GammaCGF")











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
GammaModelCGF <- function(shape, rate, iidReps = "any", ...) {
  .check_iidReps(iidReps)
  shape_fn <- validate_function_or_adaptor(shape)
  rate_fn  <- validate_function_or_adaptor(rate)
  base_cgf <- .gamma_base_cgf(iidReps = iidReps, op_name = "GammaModelCGF", ...)
  adaptCGF(
    cgf = base_cgf,
    adaptor = function(theta) c(shape_fn(theta), rate_fn(theta))
  )
}

