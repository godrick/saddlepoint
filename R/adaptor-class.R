#' @title adaptor
#'
#' @description
#' An adaptor class is an externalptr object that allows for flexible parameterization for the underlying models. It allows users to define how parameters are dynamically set or retrieved during runtime, supporting scenarios where parameters may vary or are derived from complex computations.
#'
#' There are three primary ways to create adaptors:
#' 1. **Function-based adaptors**: Use R functions to dynamically compute parameters based on input.
#' 2. **Index-based adaptors**: Subset parameters from a vector using specified indices.
#' 3. **Fixed Value adaptors**: Set parameters to known constant values.

#' @name adaptor
#' @slot .Data An externalptr object.
#' @return An adaptor object.
#'
#'
#'@noRd
setClass("adaptor", slots = list(.Data = "externalptr"))

#' Create an adaptor object
#'
#' This function constructs an adaptor object based on the provided type of argument. 
#'
#' @param indices Numeric vector of positive integers, used for creating an index-based adaptor. The indices are used to subset parameters from a vector using specified indices.
#' @param fixed_param Numeric vector of fixed parameters, used to set constant values for the parameters. 
#' @param r_func A function that dynamically compute parameters based on input vector. The function should return the required parameters based on its input.
#' 
#' @details
#' The function must be called with exactly one named argument out of `indices`, `fixed_param`, or `r_func`. 
#' For For fixed parameter values, it is recommended to use `fixed_param` argument and not `r_func` argument.
#' 
#'
#' @return An object of class "adaptor".
#' 
#' @export
#' @seealso \code{\link{BinomialModelCGF}}
#' 
#' @examples
#' \dontrun{
#'   # Example 1: Fixed and index-based parameter adaptors
#'   # Create a CGF object for a binomial distribution with a fixed parameter `n = 10`.
#'   # This configuration allows passing only the `prob` parameter in subsequent uses.
#'   binom_cgf <- BinomialModelCGF(n = adaptor(fixed_param = 10), prob = adaptor(indices = 1))
#'   # This CGF object will expect only the `prob` parameter to be passed, since `n` is already set.
#'   binom_cgf$K1(tvec = 0, parameter_vector = 0.3)
#'   
#'   # Example 3: Function adaptor
#'   # Dynamically compute probability parameter based on the input vector.
#'   r_func_cgf <- BinomialModelCGF(n = adaptor(indices = 2), prob = function(y) y[1])
#'   r_func_cgf$K1(tvec = 0, parameter_vector = c(0.3, 10))
#' }
adaptor <- function(indices = NULL, fixed_param = NULL, r_func = NULL) {
  args <- as.list(sys.call())[-1] # remove the function name from the captured call
  if (length(args) != 1) stop("Please specify exactly one of 'indices', 'fixed_param', or 'r_func'")
  if (is.null(names(args))) stop("Please specify the argument by its name")
  # if (!names(args) %in% c("indices", "fixed_param")) stop("Please specify either 'indices' or 'fixed_param'")
  
  if (names(args) == "indices"){
    if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("'indices' must be positive integers")
    return(structure(.Data = makeVectorSubsetByIndicesAdaptor(indices = indices), class = "adaptor"))
  }
  if (names(args) == "fixed_param"){
    if(!is.numeric(fixed_param)) stop(" 'fixed_param' is not defined for ", class(fixed_param)) 
    return(structure(.Data = makeSavedVectorAdaptor(fixed_parameter_values = fixed_param), class = "adaptor"))
  }
  if (names(args) == "r_func"){
    if(!is.function(r_func)) stop(" 'r_func' is not defined for ", class(r_func)) 
    if (length(formals(r_func)) != 1) stop("The function must have exactly one argument.")
    return(structure(.Data = makeAdaptorUsingRfunctions(r_function = r_func) , class = "adaptor"))
  }
}


