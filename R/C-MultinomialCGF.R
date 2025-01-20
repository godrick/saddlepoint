# R/MultinomialCGF.R


#' MultinomialCGF
#'
#' A ready-to-use CGF object for a single multinomial random vector. For a more general object that handles i.i.d. cases, use \code{\link{MultinomialModelCGF}}.
#'  
#' 
#' 
#' @details
#' **Parameter Vector Interpretation**:
#' The `parameter_vector` used when calling methods like `K1(tvec, parameter_vector)`
#' should be of the form \eqn{(N, x_1, \dots, x_d)}, where \eqn{N} is the total count,
#' and \eqn{x = x_i \dots x_d} can be any vector of non-negative entries, not all zero.
#' Therefore, \eqn{x} can be the multinomial probability vector, and can also be interpreted as a vector of odds.
#' 
#' @format An object of class \code{CGF} (an R6 class), with standard methods 
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' tvec <- rep(0,3)
#' parameter_vector <- c(10, 2, 3, 5)         # total count=10, odds=(2,3,5) or probabilities=c(0.2,0.3,0.5)
#' MultinomialCGF$K1(tvec, parameter_vector)  # MultinomialCGF$K1(tvec, c(0.2,0.3,0.5))
#'
#' @export
MultinomialCGF <- createMultinomialFamilyCGF(op_name = "MultinomialCGF")









#' Create a Parametric MultinomialCGF Object 
#'
#' @description
#' Constructs a CGF object for the multinomial distribution with a user-specified
#' parameter mapping. Optionally, you can replicate this CGF for \code{iidReps} 
#' independent and identically distributed (i.i.d.) observations.
#' 
#'
#' @param n A function of the form \code{function(theta) -> numeric}, returning the total count \eqn{n}.
#' @param prob_vec A function of the form \code{function(theta) -> numeric vector}, returning
#'   the probabilities/odds for the multinomial distribution.
#' @param block_size Either `NULL` or a positive integer specifying the size of each block
#'   for i.i.d. replication. Defaults to `NULL`.
#' @param iidReps Either `NULL` or a positive integer specifying how many i.i.d. blocks 
#'   to expect. Defaults to `NULL`.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}.
#' 
#' 
#' @return A CGF object.
#'
#'
#' @export
MultinomialModelCGF <- function(n,
                                prob_vec,
                                block_size = NULL,
                                iidReps    = NULL,
                                ...) {
  
  
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
  
  
  multinom_cgf <- createMultinomialFamilyCGF(
    op_name = "MultinomialModelCGF",
    ...
  )
  
  multinom_cgf <- adaptCGF(cgf = multinom_cgf, param_adaptor = param_adaptor)
  
  
  if (is.null(block_size) && is.null(iidReps)) return(multinom_cgf)
  if (!is.null(iidReps) && iidReps == 1) return(multinom_cgf)
  iidReplicatesCGF(cgf = multinom_cgf, iidReps = iidReps, block_size = block_size, ...)
  

}










