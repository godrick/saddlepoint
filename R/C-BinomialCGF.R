# R/C-BinomialCGF.R
# Objects: BinomialCGF, BinomialModelCGF

# Internal factory – build Binomial CGF via univariate utilities
.binomial_base_cgf <- function(iidReps, op_name, ...) {
  split_np <- function(param) {
    if (length(param) %% 2L != 0L) stop("Binomial params must be c(n_vec, p_vec) with even length.")
    M <- length(param) %/% 2L
    cbind(n = param[seq_len(M)], p = param[M + seq_len(M)])
  }

  .make_univariate_model_cgf_matrix(
    K_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      n * log(1 - p + p * exp(tvec))
    },
    K1_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      n * p * exp(tvec) / (1 - p + p * exp(tvec))
    },
    K2_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      n * p * (1 - p) * exp(tvec) / (1 - p + p * exp(tvec))^2
    },
    K3_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      tmp <- 1 - p + p * exp(tvec)
      n * (1 - p) * p * exp(tvec) * (1 - p - p * exp(tvec)) / tmp^3
    },
    K4_elem = function(tvec, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      tmp0 <- p * exp(tvec)
      tmp1 <- 1 - p
      tmp2 <- tmp1 + tmp0
      n * tmp1 * tmp0 * (4 * p * tmp0 + tmp0^2 - 4 * tmp0 + tmp1^2) / (tmp2^4)
    },
    t_hat_elem = function(x, pm) {
      n <- pm[, 1]; p <- pm[, 2]
      log((x * (1 - p)) / (p * (n - x)))
    },
    split_param_to_mat = split_np,
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
#' @examples
#' # Evaluate K1 for X ~ Binomial(10, 0.3)
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
#' `length(tvec)` must be `m * length(n(theta))`.
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
#' my_cgf <- BinomialModelCGF(n_fn, p_fn, iidReps = 1)
#' my_cgf$K1(0, c(10, 0.3))
#'
#' # non-identical example: n=c(5,10), p=c(0.2,0.7)
#' n_adapt <- function(th) th[1:2]
#' p_adapt <- function(th) th[3:4]
#' my_cgf2 <- BinomialModelCGF(n_adapt, p_adapt, iidReps = "any")
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
    if (length(n_vals) != length(p_vals)) stop("n(theta) and prob(theta) must have the same length.")
    c(n_vals, p_vals)
  }

  adaptCGF(base_cgf, adaptor_fun)
}

