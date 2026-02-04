# R/SubunitaryMultinomialCGF.R
# Objects: SubunitaryMultinomialCGF, SubunitaryMultinomialModelCGF
#
# This file defines a "subunitary" multinomial CGF: a multinomial family where
# one extra category is known to have zero count.
#
# With the current multinomial implementation (MultinomialFamilyCGF +
# createMultinomialFamilyCGF), replication across i.i.d. blocks is handled by
# createMultinomialFamilyCGF() via iidReplicatesCGF().
#
# Therefore, the functions below implement ONLY the single-block formulas:
# they assume length(tvec) == d, where d = length(parameter_vector) - 1.
#
# The subunitary CGF arises by starting from a (d+1)-category multinomial and
# setting the (d+1)th CGF argument to -Inf, which effectively multiplies by the
# indicator 1{Y_{d+1}=0}. This gives:
#
#   K_sub(t_1,...,t_d) = N * log( sum_{i=1}^d pi_i * exp(t_i) )
#
# where parameter_vector = c(N, pi_1,...,pi_d) and sum(pi) <= 1.
#
# Compared to the standard multinomial family (which uses probabilities
# normalized as p_i = pi_i / sum(pi)), K_sub differs only by the additive
# constant N * log(sum(pi)). That constant affect likelihood derivatives
# with respect to the model parameters, but it does NOT change K1/K2/K3/K4
# as derivatives with respect to t.
#
# We therefore implement SubunitaryMultinomial* by overriding only K() and
# inheriting all higher-order t-derivative logic from MultinomialFamilyCGF.

# -----------------------------------------------------------------------------
# Single-block K override
# -----------------------------------------------------------------------------

#' @noRd
.subunitaryMultinomial_K_default <- function(tvec, parameter_vector) {

  d <- length(parameter_vector) - 1L
  if (d < 1L) {
    stop("SubunitaryMultinomial: 'parameter_vector' must have length >= 2 (c(N, pi[1:d])).")
  }
  if (length(tvec) != d) {
    stop(
      "SubunitaryMultinomial::K: length(tvec) must equal d = length(parameter_vector) - 1. ",
      "Got length(tvec) = ", length(tvec), ", d = ", d, "."
    )
  }

  N_val  <- parameter_vector[1]
  pi_val <- parameter_vector[-1]
  pi_sum <- sum(pi_val)


  # we use exp(t) - 1 rather than expm1(t) # RTMB-friendly!!!
  zm1  <- exp(tvec) - 1
  frac <- sum(pi_val * zm1) / pi_sum

  # K_sub(t) = N*log(pi_sum) + N*log1p( sum((pi/pi_sum)*(exp(t)-1)) )
  #          = N*log( sum_i pi_i exp(t_i) )
  N_val * (log(pi_sum) + log1p(frac))
}



# -----------------------------------------------------------------------------
# Ready-to-use CGF object
# -----------------------------------------------------------------------------

#' SubunitaryMultinomialCGF
#'
#' A ready-to-use CGF object for the subunitary multinomial model.
#'
#' Parameter vector: c(N, pi_1, ..., pi_d), where pi_i are the probabilities
#' of the 'active' categories and sum(pi) <= 1. The omitted category has
#' probability 1 - sum(pi) and is assumed to have count 0.
#'
#' Technically, this is the (d+1)-category multinomial CGF evaluated at
#' \eqn{t_{d+1} = -Inf}, which introduces a multiplicative factor for the event
#' \eqn{\{Y_{d+1}=0\}}. This is why the K() differs from the standard multinomial
#' family by an additive constant N*log(sum(pi)).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' d <- 3
#' N <- 10
#' pi_vec <- c(0.2, 0.3, 0.1)  # sum < 1 (subunitary)
#' tt <- rnorm(d)
#'
#' param <- c(N, pi_vec)
#'
#' # K1 (and all higher t-derivatives) match the standard multinomial family
#' SubunitaryMultinomialCGF$K1(tt, param)
#' MultinomialCGF$K1(tt, param)
#'
#'
#' # K differs by N*log(sum(pi_vec)) per iid block:
#' SubunitaryMultinomialCGF$K(tt, param) - MultinomialCGF$K(tt, param)
#' # should be ~ N * log(sum(pi_vec))
#'
#' # With multiple iid blocks,
#' # the difference scales with the number of blocks:
#' t2 <- c(tt, tt)
#' SubunitaryMultinomialCGF$K(t2, param) - MultinomialCGF$K(t2, param)
#' # should be ~ 2*N*log(sum(pi_vec))
#' }

#'
#' @format An object of class CGF (R6), with methods K, K1, K2, K3operator, K4operator, etc.
#' @export
SubunitaryMultinomialCGF <- createMultinomialFamilyCGF(
  op_name = "SubunitaryMultinomialCGF",
  K = .subunitaryMultinomial_K_default
)



# -----------------------------------------------------------------------------
# Parametric (theta-mapped) constructor
# -----------------------------------------------------------------------------

#' Create a Parametric Subunitary Multinomial CGF Object
#'
#' @description
#' Creates a CGF object for the subunitary multinomial distribution with
#' total count N = n(theta) and active-category probabilities pi = prob_vec(theta),
#' where sum(pi) <= 1.
#'
#'
#' @details
#' Let \eqn{Y \sim \mathrm{Multinomial}(N,\pi_1,\dots,\pi_d,\pi_{d+1})} with
#' \eqn{\sum_{i=1}^{d+1} \pi_i = 1}. Conditioning on \eqn{\{Y_{d+1}=0\}}, we have
#' \eqn{W = (Y \mid Y_{d+1}=0) \sim \mathrm{Multinomial}(N,p_1,\dots,p_d) }, where
#' \deqn{
#'   p_i = \frac{\pi_i}{\sum_{j=1}^d \pi_j}, \quad i=1,\dots,d.
#' }
#'
#'
#' **The CGF:**
#'
#' From the perspective of restricting to \eqn{\{Y_{d+1}=0\}}, one effectively sets
#' \eqn{t_{d+1}=-\infty} in the multinomial CGF. This yields the \emph{subunitary} CGF:
#' \deqn{
#'   K_{\mathrm{subunitary}}(t_1,\dots,t_d)
#'   = \log\bigl\{\Pr(Y_{d+1}=0)\bigr\}
#'   \;+\;
#'   K_{W}(t_1,\dots,t_d)\,,
#' }
#' where \eqn{\Pr(Y_{d+1}=0) = N \log\left(\!\sum_{i=1}^d \pi_i\right)} and
#' \eqn{K_{W}} is the usual multinomial CGF for \eqn{\mathrm{Multinomial}(N,p_1,\dots,p_d)}.
#'
#' In many applications (such as certain capture-recapture models), some categories
#' are effectively impossible based on observed data, so they can be merged into
#' a single "impossible" category with zero count. The probabilities of the
#' \eqn{d} active categories then sum to less than one.
#'
#'
#' @param n A function (or adaptor) mapping theta to the multinomial total count N.
#' @param prob_vec A function (or adaptor) mapping theta to the active-category
#'   probability vector pi (length d). These entries should be non-negative and
#'   sum(pi) should be <= 1.
#' @param iidReps Either "any" or a positive integer specifying how many i.i.d.
#'   blocks are expected.
#' @param ... Additional named arguments passed to CGF creating functions.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' d <- 3
#' N <- 10
#' pi_vec <- c(0.2, 0.3, 0.1)  # sum < 1
#' theta <- c(N, pi_vec)
#'
#' # Build a theta-mapped CGF: theta -> (N, pi_1,...,pi_d)
#' cgf <- SubunitaryMultinomialModelCGF(
#'   n        = function(th) th[1],
#'   prob_vec = function(th) th[2:(d + 1)],
#'   iidReps  = 2
#' )
#'
#' # With iidReps = 2, tvec must be length 2*d:
#' t_block <- rnorm(d)
#' tvec <- c(t_block, t_block)
#'
#' # Evaluate K and K1 at theta:
#' cgf$K(tvec, theta)
#' cgf$K1(tvec, theta)
#' }
#'
#' @return A CGF object.
#' @export
SubunitaryMultinomialModelCGF <- function(n,
                                         prob_vec,
                                         iidReps = "any",
                                         ...) {

  .check_iidReps(iidReps)

  sub_cgf <- createMultinomialFamilyCGF(
    iidReps = iidReps,
    op_name = "SubunitaryMultinomialModelCGF",
    K = .subunitaryMultinomial_K_default,
    ...
  )

  n_fn <- validate_function_or_adaptor(n)
  prob_vec_fn <- validate_function_or_adaptor(prob_vec)

  param_adaptor <- function(theta) c(n_fn(theta), prob_vec_fn(theta))

  adaptCGF(cgf = sub_cgf, adaptor = param_adaptor)
}
