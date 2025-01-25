# R/MultinomialCGF.R


#' MultinomialCGF
#' 
#' A ready-to-use CGF object for handling multiple i.i.d. multinomial random vectors. 
#' The number of multinomial random vectors is determined by the length of `tvec` and the length of the multinomial probability vector.
#' For user-specified parameter mappings, refer to \code{\link{MultinomialModelCGF}}.
#' 
#' 
#' @details
#' **Parameter Vector Interpretation**:
#' The `parameter_vector` used when calling methods like \code{K1(tvec, parameter_vector)}
#' should be of the form \eqn{(n, x_1, \dots, x_d)}, where:
#' \itemize{
#'   \item \eqn{n} is the total count (a positive integer),
#'   \item \eqn{x = (x_1, \dots, x_d)} is a vector of non-negative entries, not all zero.
#' }
#' 
#' The vector \eqn{x} can be interpreted in two ways:
#' \itemize{
#'   \item **Odds**: \eqn{x_i} are odds values. Probabilities are derived as \eqn{p_i = x_i / \sum_{j=1}^{d} x_j}.
#'   \item **Probabilities**: If \eqn{x} sums to 1, it is treated directly as the probability vector \eqn{p}.
#' }
#' 
#' **IID Assumption**:
#' All multinomial random vectors processed by `MultinomialCGF` are assumed to be i.i.d. 
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
#' Constructs a CGF object for the multinomial distribution based on user-specified
#' functions or adaptors that map a parameter vector \code{theta} to the multinomial's total count \eqn{n}
#' and probability vector \eqn{p}. Optionally, you can replicate this CGF for \code{iidReps} 
#' i.i.d. multinomial observations, thus imposing a length restriction on \code{tvec}.
#' 
#' 
#' @details
#' **User-Specified Parameter Mapping**:
#' \itemize{
#'   \item \code{n}: An adaptor or a function of the form \code{function(theta) -> numeric}. 
#'         Must return a single numeric value, interpreted as the total count \eqn{n} for the multinomial distribution.
#'         The returned value should be non-negative.
#'   \item \code{prob_vec}: An adaptor or a function of the form \code{function(theta) -> numeric vector}. 
#'         Must return a vector of non-negative entries, interpreted as probabilities or odds for the multinomial distribution.
#'         If \eqn{prob_vec} sums to 1, it is treated as a probability vector; otherwise, it is normalized into probabilities.
#' }
#' 
#' **I.I.D. Replicates**:
#' By setting \code{iidReps} to a positive integer \eqn{m}, you declare that the input vector \eqn{tvec} will be split
#' into \eqn{m} blocks of equal size, each corresponding to one i.i.d. multinomial sample. 
#' If \code{iidReps} is \code{"any"}, no length restriction is enforced on \eqn{tvec}, allowing flexible usage where the dimension 
#' of \eqn{tvec} can vary without implying i.i.d. blocks.
#' 
#' @param n An adaptor or function defining the total count \eqn{n(\theta)}. 
#'   It must accept a parameter vector \eqn{\theta} and return a single numeric value.
#' @param prob_vec An adaptor or function defining the probability (or odds) vector \eqn{p(\theta)}. 
#'   It must accept a parameter vector \eqn{\theta} and return a numeric vector.
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, such as method overrides or operator definitions.
#' 
#' @return A `CGF` object (an R6 class) specialized for the multinomial distribution with user-defined mappings.
#'   The returned object supports the usual CGF methods (\code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.).
#'   
#' @examples
#' # Suppose we want n = 10, and probabilities = c(0.2, 0.3, 0.5).
#' n_adaptor <- adaptor(fixed_param = 10)
#' p_adaptor <- adaptor(fixed_param = c(0.2, 0.3, 0.5))
#' 
#' my_cgf <- MultinomialModelCGF(n = n_adaptor, prob_vec = p_adaptor, iidReps = 2)
#' 
#' # The resulting CGF object expects 2 i.i.d. multinomial blocks,
#' # each of dimension 3. So tvec must be of length 6.
#' tvec <- rep(0, 6)
#' param <- c(0.2, 0.3, 0.5)  # or any unused parameter vector
#' 
#' # We can now compute, e.g., the gradient:
#' my_cgf$K1(tvec, param)
#'
#'
#' @export
MultinomialModelCGF <- function(n,
                                prob_vec,
                                iidReps = "any",
                                ...) {
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  
  multinom_cgf <- createMultinomialFamilyCGF(
    op_name = "MultinomialModelCGF",
    ...
  )
  
  n_fn <- validate_function_or_adaptor(n)
  prob_vec_fn <- validate_function_or_adaptor(prob_vec)
  # adapt the user param into c(n, p1, p2, ..., p_d) 
  param_adaptor <- function(theta) c(n_fn(theta), prob_vec_fn(theta))
  
  adaptCGF(cgf = multinom_cgf, param_adaptor = param_adaptor)
}










