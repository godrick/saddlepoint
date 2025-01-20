# R/SubunitaryMultinomialCGF.R



#' @title SubunitaryMultinomialModelCGF
#'
#' @description
#' Generates a CGF object derived from a modified Multinomial distribution. This modification accommodates instances where certain outcomes are explicitly zero, and consequently, the probabilities associated with the remaining non-zero outcomes sum to a value less than one.
#'
#' @details
#' Given \eqn{x_1,...,x_k} as the non-zero outcomes of a Multinomial distribution that sum to \eqn{N},
#' with each outcome having an associated probability \eqn{p_i} where \eqn{ \sum_{i=1}^{k} p_i < 1},
#' this CGF is designed to integrate the occurrence of a known zero outcome, \eqn{x_{k+1}}, thereby ensuring \eqn{ \sum_{i=1}^{k+1} p_i = 1}.
#'
#' This approach allows for the modeling and analysis of Multinomial distributions in situations where some categories are known not to be observed and the probabilities of observing the other categories are adjusted accordingly.
#'
#'
#' @param n An \code{adaptor} or a function of the form \code{function(theta) -> numeric}, returning the number of trials in the multinomial distribution.
#' @param prob_vec An \code{adaptor} or a function of the form \code{function(theta) -> numeric vector}, returning
#'   the multinomial probabilities. The sum of these probabilities should be strictly less than one to validate the occurrence of the known zero outcome.
#' @param block_size Either `NULL` or a positive integer specifying the size of each block
#'   for i.i.d. replication. Defaults to `NULL`.
#' @param iidReps Either `NULL` or a positive integer specifying how many i.i.d. blocks 
#'   to expect. Defaults to `NULL`.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}.
#'
#'
#' @return A \eqn{CGF} object.
#' @seealso \code{\link{MultinomialModelCGF}}
#'
#' @export
SubunitaryMultinomialModelCGF <- function(n, 
                                          prob_vec,
                                          block_size = NULL,
                                          iidReps    = NULL,
                                          ...){
  # Only one of (iidReps, block_size) can be set:
  if (!is.null(iidReps) && !is.null(block_size)) stop("Please specify only one of 'iidReps' or 'block_size', not both.")
  
  # if iidReps is set:
  if (!is.null(iidReps)) {
    if (!is.numeric(iidReps) || length(iidReps) != 1 ||
        iidReps < 1 || iidReps != as.integer(iidReps)) {
      stop("'iidReps' must be NULL or a positive integer.")
    }
  }
  
  # if block_size is set:
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1 ||
        block_size < 1 || block_size != as.integer(block_size)) {
      stop("'block_size' must be NULL or a positive integer.")
    }
  }
  
  
  n_fn <- validate_function_or_adaptor(n)
  prob_vec_fn <- validate_function_or_adaptor(prob_vec)
  
  # adapt the user param into c(n, p1, p2, ..., p_d) 
  param_adaptor <- function(theta) c(n_fn(theta), prob_vec_fn(theta))
  
  
  K_func_default <- function(tvec, param) {
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
    n_val*log(prob_sum) + n_val*log1p( sum(adjusted_prob_vec*zm1) )
  }
  
  multinom_cgf <- createMultinomialFamilyCGF(
    K_func = K_func_default,
    op_name = "SubunitaryMultinomialModelCGF",
    ...
  )
  
  multinom_cgf <- adaptCGF(cgf = multinom_cgf, param_adaptor = param_adaptor)
  
  
  if (is.null(block_size) && is.null(iidReps)) return(multinom_cgf)
  if (!is.null(iidReps) && iidReps == 1) return(multinom_cgf)
  iidReplicatesCGF(cgf = multinom_cgf, iidReps = iidReps, block_size = block_size, ...)
}