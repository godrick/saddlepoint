# R/SubunitaryMultinomialCGF.R


.subunitaryMultinomial_K_func_default <- function(tvec, param) {
  # The CGF for a modified Multinomial distribution.
  # It is useful in scenarios where a multinomial category are known to be 0.
  # 
  # For a random vector Y = (Y_1,...,Y_{d+1}),
  # prob_vec represents the probability vector for the first 'd' categories in a Multinomial distribution.
  # These probabilities can have a sum less than 1, indicating the presence of an additional zero-outcome category.
  # 
  # tvec is the vector of CGF arguments corresponding to each category, with the last category (d+1) assumed to be zero.
  # N is the parameter of the Multinomial distribution, representing the number of trials.
  # 
  # The function computes the CGF by adjusting the probability vector to account for the zero-outcome category.
  # The adjusted probabilities are calculated as p_i = prob_vec[i] / sum(prob_vec), for i = 1 to d.
  # K = N*log(sum_i prob_vec_i) + N*log(sum_i p_i e^{t_i})
  # 
  # zm1, the vector whose entries are zm1[i] = exp(tvec[i])-1
  
  n_val    <- param[1]
  prob_vec <- param[-1]
  
  zm1 <- exp(tvec) - 1
  prob_sum <- sum(prob_vec)
  adjusted_prob_vec = prob_vec / prob_sum
  # n_val*log(prob_sum) + n_val*log1p( sum(adjusted_prob_vec*zm1) )
  
  
  
  d <- length(adjusted_prob_vec) 
  nblocks <- length(tvec) / d
  b_zm1 <- matrix(zm1, nrow = d, ncol = nblocks)
  block_results <- apply(b_zm1, 2, function(zm1_block) {
    n_val*log(prob_sum) + n_val*log1p( sum(adjusted_prob_vec*zm1_block) )
  })
  sum(block_results)
}




#' @title Subunitary Multinomial CGF
#'
#' @description
#' Constructs the CGF object for a multinomial distribution 
#' in which one of the categories is \emph{known} to have zero counts, leaving 
#' the remaining \eqn{d} categories with a total probability less than one.
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
#'
#'
#' @param n An `adaptor` or a function of the form \code{function(theta) -> numeric}, returning the number of trials in the multinomial distribution.
#' @param prob_vec An `adaptor` or a function of the form \code{function(theta) -> numeric}, 
#'   returning the multinomial probabilities for the \eqn{d} active categories. 
#'   The sum of these probabilities should be strictly less than one to account for 
#'   the excluded category.
# #' @param iidReps Either \code{"any"} or a positive integer specifying how many
# #'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, such as method overrides or operator definitions.
#'
#' @return A `CGF` object.
#' @seealso \code{\link{MultinomialModelCGF}}
#'
#' @export
SubunitaryMultinomialModelCGF <- function(n, 
                                          prob_vec,
                                          ...){
  
  
  multinom_cgf <- createMultinomialFamilyCGF(
    K_func = .subunitaryMultinomial_K_func_default,
    op_name = "SubunitaryMultinomialModelCGF",
    ...
  )
  
  n_fn <- validate_function_or_adaptor(n)
  prob_vec_fn <- validate_function_or_adaptor(prob_vec)
  
  param_adaptor <- function(theta) c(n_fn(theta), prob_vec_fn(theta))
  
  adaptCGF(cgf = multinom_cgf, param_adaptor = param_adaptor)
  
}