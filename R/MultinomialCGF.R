# R/MultinomialCGF.R


#' MultinomialCGF
#'
#' A ready-to-use CGF object for the multinomial distribution.
#' 
#' @details
#' **Parameter Vector Interpretation**:
#' The `parameter_vector` used when calling methods like `K(tvec, parameter_vector)`
#' should be of the form \eqn{(N, x_1, \dots, x_d)}, where \eqn{N} is the total count,
#' and \eqn{x = x_i \dots x_d} can be any vector of non-zero numbers not all zero.
#' Therefore, \eqn{x} can be the multinomial probability vector, and can also be interpreted as a vector of odds.
#' 
#' @format An object of class \code{CGF} (an R6 class), with standard methods 
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, \code{K4operator}, etc.
#'
#' @examples
#' tvec <- c(0.1, -0.2, 0.05)
#' param <- c(10, 2, 3, 5) # total count=10, odds=(2,3,5)
#' MultinomialCGF$K(tvec, param)
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
#' @param N A function of the form \code{function(theta) -> numeric}, returning the total count \eqn{N}.
#' @param prob_vec A function of the form \code{function(theta) -> numeric vector}, returning
#'   the probabilities/odds for the multinomial distribution.
#' @param iidReps Integer. If \code{iidReps > 1}, we replicate the resulting CGF for 
#'   multiple i.i.d. copies of the multinomial variables. Defaults to \code{1}.
#' @param ... Additional arguments passed to \code{\link{createCGF}} or possibly
#'   to \code{\link{iidReplicatesCGF}} if \code{iidReps > 1}.
#' 
#' 
#' @return A CGF object.
#'
#'
#' @export
MultinomialModelCGF <- function(N,
                                prob_vec,
                                iidReps = 1,
                                ...) {
  
  # # quick checks
  # if (!is.function(N)) stop("`N` must be a valid function that returns the total count parameter of the multinomial distribution.")
  # if (!is.function(prob_vec)) stop("`prob_vec` must be a valid function that returns the probability vector of the multinomial distribution.")
  # if (length(formals(N)) != 1) stop("`N` must be a function that accepts exactly one argument.")
  # if (length(formals(prob_vec)) != 1) stop("`prob_vec` must be a function that accepts exactly one argument.")
  
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1) stop("'iidReps' must be a positive integer.")
  N_fn <- validate_function_or_adaptor(N)
  prob_vec_fn <- validate_function_or_adaptor(prob_vec)
  
  # adapt the user param into c(N, p1, p2, ..., p_d) 
  param_adaptor <- function(theta) c(N_fn(theta), prob_vec_fn(theta))
  
  
  multinom_cgf <- createMultinomialFamilyCGF(
    param_adaptor = param_adaptor,
    op_name = "MultinomialModelCGF",
    ...
  )
  
  if (iidReps == 1) return(multinom_cgf)
  
  iidReplicatesCGF(cgf = multinom_cgf, iidReps = iidReps, ...)
}
