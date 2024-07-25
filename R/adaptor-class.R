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

#' adaptor constructor
#'
#' This function constructs an object of class "adaptor" using a valid externalptr object.
#'
#' @param ptr An externalptr object. It should not be an empty externalptr object.
#'
#' @return An object of class "adaptor".
#' @examples
#' \dontrun{
#'   # assuming ptr is a valid externalptr object
#'   adaptor <- adaptor(ptr)
#' }
#' @noRd
adaptor <- function(ptr) {
  if (!is(ptr, "externalptr")) stop("ptr must be of type externalptr")
  if (isTRUE( identical(ptr, new("externalptr"))) ) stop("ptr is an empty externalptr object")
  structure(.Data = ptr, class = "adaptor")
}


#' @title Create fixed parameter adaptor
#'
#' @description
#' Instantiates an adaptor with fixed parameter values. Useful for models where some parameters do not change.
#'
#'
#' @param fixed_param Numeric vector of fixed parameters.
#'
#' @return An adaptor object set with fixed parameters.
#'
#' @examples
#' \dontrun{
#' # This is especially useful for setting known values of parameters that will not change during estimation.
#' # For instance, one might set the `n` parameter to `1` for a Bernoulli distribution (which is a special case of the binomial distribution).
#' }
#'
#' @export
#' @seealso \code{\link{indices.adaptor}}
fixed.parameter.adaptor <- function(fixed_param){
  if(!is.numeric(fixed_param)) stop(" 'fixed_param' is not defined for ", class(fixed_param)) 
  adaptor(makeSavedVectorAdaptor(fixed_parameter_values = fixed_param))
}



#' @title Create index-based adaptor
#'
#' @description
#' Creates an adaptor that subsets parameters from a vector using specified indices.
#'
#' @param indices Numeric vector of positive integers for subsetting parameters.
#'
#' @return An adaptor object that subsets parameters by indices.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' @seealso \code{\link{fixed.parameter.adaptor}}
indices.adaptor <- function(indices){
  if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("indices must be positive integers")
  adaptor(makeVectorSubsetByIndicesAdaptor(indices = indices))
}
