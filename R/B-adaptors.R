# R/adaptors.R
# Objects: adaptor class and its constructor, validate_function_or_adaptor








#' @title adaptor
#'
#' @description
#' An S4 class used to encapsulate how a parameter is extracted or computed
#' from a parameter vector \code{theta}. This can be:
#' \itemize{
#'   \item A subset of \code{theta} by index,
#'   \item A fixed numeric value,
#'   \item A user-defined function.
#' }
#'
#' Internally, the `func` slot is a plain R function of the form
#' \code{function(param_vector) -> numeric}.
#'
#' @slot func A function(\code{param_vector}) -> numeric, implementing the logic.
#' @noRd
setClass("adaptor", slots = list(func = "function"))





#' Create an `adaptor` object
#'
#' This function constructs an \code{adaptor} object from exactly one of:
#' \enumerate{
#'   \item \code{indices}: A numeric vector of positive integer indices, used to subset
#'     \code{theta}.
#'   \item \code{fixed_param}: A numeric value (or vector) to be returned regardless of \code{theta}.
#'   \item \code{r_func}: A custom R function of the form
#'     \code{function(theta) -> numeric}.
#' }
#'
#' @param indices A numeric vector of positive integers.
#' @param fixed_param A numeric vector or scalar to return.
#' @param r_func A function(\code{theta}) -> numeric.
#'
#' @details
#' Exactly **one** of \code{indices}, \code{fixed_param}, or \code{r_func} must be non-\code{NULL}.
#'
#'
#' @return An object of class \code{adaptor}.
#'
#'
#' @export
#' @examples
#' \dontrun{
#'   # Example 1: Fixed and index-based parameter adaptors
#'   # Create a CGF object for a binomial distribution with a fixed parameter `n = 10`.
#'   # This configuration allows passing only the `prob` parameter in subsequent uses.
#'   binom_cgf <- BinomialModelCGF(n = adaptor(fixed_param = 10), prob = adaptor(indices = 1))
#'   # The object binom_cgf will expect only the `prob` parameter to be passed, since `n` is already set.
#'   binom_cgf$K1(tvec = 0, parameter_vector = 0.3)
#'   
#'   # Example 2: Function adaptor
#'   # Dynamically compute probability parameter based on the input vector.
#'   cgf <- BinomialModelCGF(n = adaptor(indices = 2), prob = function(y) y[1])
#'   cgf$K1(tvec = 0, parameter_vector = c(0.3, 10))
#' }
adaptor <- function(indices = NULL, fixed_param = NULL, r_func = NULL) {
  args <- as.list(sys.call())[-1] # remove the function name from the captured call
  if (length(args) != 1) stop("Please specify exactly one of 'indices', 'fixed_param', or 'r_func'")
  if (is.null(names(args))) stop("Please specify the argument by its name, e.g. `indices = c(...)`.")


  if (names(args) == "indices") {
    if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("'indices' must be positive integers")
    func_ <- function(param_vector) {
      param_vector[indices]
    }
    return(methods::new("adaptor", func = func_))
  }

  if (names(args) == "fixed_param"){
    if(!is.numeric(fixed_param)) stop("'fixed_param' must be numeric.")
    func_ <- function(param_vector) { fixed_param + 0*param_vector[1] }
    return(methods::new("adaptor", func = func_))
  }

  if (names(args) == "r_func"){
    if(!is.function(r_func))  stop("'r_func' must be a function of form function(theta)->numeric.")
    if (length(formals(r_func)) != 1)  stop("'r_func' must have exactly one argument (e.g. function(theta)).")
    return(methods::new("adaptor", func = r_func))
  }

  stop("Internal error: no valid argument chosen. Check your usage of `adaptor()`.")
}






# #' Evaluate an \code{adaptor} or an R function with a parameter vector
# #'
# #' @param obj An object of class \code{adaptor} or an R function.
# #' @param param_vector A numeric vector of parameters.
# #' @return The numeric result of calling \code{ad@.Data(param_vector)}.
# #' @noRd
# check_and_evaluate_adaptor <- function(obj, param_vector) {
#   if (is.function(obj)) {
#     if (length(formals(obj)) != 1) stop("The function must accept exactly one argument.")
#     return(obj(param_vector))
#   } else if (inherits(obj, "adaptor")) {
#     fn <- obj@.Data
#     return(fn(param_vector))
#   } else {
#     stop("Not a valid 'adaptor' object.")
#   }
# }


#' Validate an `adaptor` or an R function, returning a function
#' 
#' This helper function ensures that the provided object is either:
#' \enumerate{
#'   \item A user-defined function(\code{theta}) -> numeric,
#'   \item An `adaptor` object whose `func` slot is such a function.
#' }
#'
#' @param obj Either a function(\code{theta}) -> numeric or an object of class `adaptor`.
#' @return A function of the form \code{function(theta) -> numeric}.
#' @noRd
validate_function_or_adaptor <- function(obj) {
  if (!is.function(obj)) {
    if (inherits(obj, "adaptor")) {
      return(obj@func)
    }
    stop("Provided argument must be either a function or an object of type 'adaptor'.")
  }
  if (length(formals(obj)) != 1) stop("The function provided must have exactly one argument (e.g. function(theta)).")
  obj
}






