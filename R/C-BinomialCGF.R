# R/C-BinomialCGF.R
# Objects: BinomialCGF, BinomialModelCGF

# Internal factory – build Binomial CGF via univariate utilities
.binomial_base_cgf <- function(iidReps, op_name, ...) {
  split_np <- function(param) {
    ln <- length(param) / 2
    if (ln != as.integer(ln) || ln < 1) stop("Binomial parameters must be concatenated as c(n[1:L], p[1:L]).")
    idx <- seq_len(ln)
    cbind(n = param[idx], p = param[ln + idx])
  }

  .make_univariate_model_cgf_matrix(
    K_elem = function(tvec, pm) { pm[,1] * log(1 - pm[,2] + pm[,2]*exp(tvec))},
    K1_elem = function(tvec, pm) { pm[,1] * pm[,2] * exp(tvec) / (1 - pm[,2] + pm[,2]*exp(tvec)) },
    # NOTE (AD / Hessian safety):
    # we avoid "x^k" on quantities that may be negative under AD (higher-order derivatives can produce NaN in RTMB).
    K2_elem = function(tvec, pm) {
      denom <- 1 - pm[,2] + pm[,2] * exp(tvec)
      denom2 <- denom * denom
      pm[,1] * pm[,2] * (1 - pm[,2]) * exp(tvec) / denom2
    },

    # K4_elem = function(tvec, pm) {
    #   n <- pm[, 1]; p <- pm[, 2]
    #   tmp0 <- p * exp(tvec)
    #   tmp1 <- 1 - p
    #   tmp2 <- tmp1 + tmp0
    #   n * tmp1 * tmp0 * (4 * p * tmp0 + tmp0^2 - 4 * tmp0 + tmp1^2) / (tmp2^4)
    # },

    K3_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      tmp <- 1 - p + p * exp(tvec)
      tmp2 <- tmp * tmp
      tmp3 <- tmp2 * tmp
      n * (1 - p) * p * exp(tvec) * (1 - p - p * exp(tvec)) / tmp3
    },
    K4_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      tmp0 <- p * exp(tvec)
      tmp1 <- 1 - p
      tmp2 <- tmp1 + tmp0
      tmp0sq <- tmp0 * tmp0
      tmp1sq <- tmp1 * tmp1
      tmp2sq <- tmp2 * tmp2
      tmp2four <- tmp2sq * tmp2sq
      n * tmp1 * tmp0 * (4 * p * tmp0 + tmp0sq - 4 * tmp0 + tmp1sq) / tmp2four
    },

    t_hat_elem = function(x, pm) {
      log(x * (1 - pm[,2])) - log(pm[,2] * (pm[,1] - x))
    },
    split_param_to_mat = split_np,
    simulate_func = function(iidReps, parameter_vector, ...) {
      pm <- split_np(parameter_vector)
      n <- as.numeric(pm[, 1])
      p <- as.numeric(pm[, 2])

      if (any(!is.finite(n)) || any(n < 0)) {
        stop("BinomialCGF$rsim: 'n' must be finite and >= 0.")
      }
      # For simulation, n must be integer-ish
      if (any(abs(n - round(n)) > 1e-8)) stop("BinomialCGF$rsim: 'n' must be an integer to simulate.")
      n <- as.integer(round(n))
      if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("BinomialCGF$rsim: 'p' must be finite and in [0, 1].")

      d <- length(n)
      out <- stats::rbinom(
        n    = d * iidReps,
        size = rep.int(n, times = iidReps),
        prob = rep.int(p, times = iidReps)
      )
      matrix(out, nrow = d, ncol = iidReps)
    },

    iidReps = iidReps,
    op_name = op_name,
    ...
  )
}

#' Binomial CGF Object
#'
#' Ready-to-use CGF for Binomial(n, p). Accepts packed params `c(n_vec, p_vec)`.
#' If `length(tvec)` is a multiple of `length(n_vec)`, evaluation proceeds (iidReps = "any").
#'
#' By default, \code{BinomialCGF} supports vectorized evaluation for i.i.d. replicates.
#'
#' @format An object of class \code{CGF} (R6), with standard methods:
#'   \code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.
#'
#' @examples
#' # Expected value of X ~ Binomial(10, 0.3)
#' BinomialCGF$K1(0, c(10, 0.3))
#'
#' @export
BinomialCGF <- .binomial_base_cgf(
  iidReps = "any",
  op_name = "BinomialCGF"
)

#' Create a Parametric Binomial CGF Object
#'
#' @description
#' Binomial with parameters `n(theta)` and `prob(theta)`. If `iidReps = "any"` (default),
#' `length(tvec)` must be a multiple of `length(n(theta))`. If `iidReps = m`, then
#' `length(tvec)` must be `m * length(n(theta))`. It supports i.i.d. replication (via \code{iidReps}) and non-identical usage by
#' causing `n(theta)` and `prob(theta)` to return vectors.
#'
#' @param n A function or `adaptor` mapping `theta` -> scalar or vector of trial counts.
#' @param prob A function or `adaptor` mapping `theta` -> scalar or vector of probabilities.
#' @param iidReps Either `"any"` or a positive integer. Default `"any"`.
#' @param ... Passed through to the CGF creation (e.g., to override operators).
#' @return A `CGF` object.
#'
#' @examples
#' # i.i.d. example: n=10, p=0.3
#' n_fn <- function(th) th[1]
#' p_fn <- function(th) th[2]
#' my_cgf <- BinomialModelCGF(n_fn, p_fn, iidReps = 1) # if iidReps="any", length(tvec) determines the number of i.i.d. replicates
#' my_cgf$K1(0, c(10, 0.3))
#'
#' # non-identical example: n=c(5,10), p=c(0.2,0.7)
#' # param ==> c(5,10, 0.2, 0.7)
#' n_adapt <- function(th) th[1:2] # OR n_adapt <- adaptor(indices = 1:2)
#' p_adapt <- function(th) th[3:4] # OR p_adapt <- adaptor(indices = 3:4)
#' my_cgf2 <- BinomialModelCGF(n_adapt, p_adapt, iidReps="any") # default iidReps="any", if iidReps=m, then m-i.i.d. blocks are expected => length(tvec) must be a multiple of m
#' # e.g. tvec= c(0,0) => length=2 => we have 2 binomial distributions: (5,0.2) and (10,0.7)
#' my_cgf2$K1(c(0,0), c(5,10,0.2,0.7))
#'
#' @export
BinomialModelCGF <- function(n, prob, iidReps = "any", ...) {
  .check_iidReps(iidReps)

  n_fn <- validate_function_or_adaptor(n)
  p_fn <- validate_function_or_adaptor(prob)

  base_cgf <- .binomial_base_cgf(iidReps = iidReps, op_name = "BinomialModelCGF", ...)

  adaptor_fun <- function(theta) {
    n_vals <- n_fn(theta)
    p_vals <- p_fn(theta)
    if (length(n_vals) != length(p_vals)) stop( sprintf("Length mismatch: length(n)=%d, length(p)=%d.", length(n_vals), length(p_vals)) )
    c(n_vals, p_vals)
  }

  adaptCGF(base_cgf, adaptor_fun)
}

