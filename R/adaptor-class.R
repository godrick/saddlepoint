#' @title adaptor
#'
#' @description
#' The adaptor class is an externalptr object that points to related C++ objects.
#'
#' @slot .Data An externalptr object.
#' @return An adaptor object.
#'
#'
#'@noRd
setClass("adaptor", slots = list(.Data = "externalptr"))

#' adaptor Constructor
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


#' @title Create an adaptor using fixed values
#'
#' @description
#' Creates an adaptor using a specified numeric vector of fixed parameter values.
#' This is especially useful for setting known values of parameters that will not change during estimation.
#' For instance, one might set the `n` parameter to `1` for a Bernoulli distribution (which is a special case of the binomial distribution).
#'
#'
#' @param fixed_param A numeric vector indicating the fixed parameter values.
#'
#' @return "adaptor" object
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' @seealso \code{\link{indices_adaptor}}
fixed_parameter_adaptor <- function(fixed_param){
  if(!is.numeric(fixed_param)) stop(" 'fixed_param' is not defined for ", class(fixed_param)) 
  adaptor(makeSavedVectorAdaptor(fixed_parameter_values = fixed_param))
}



#' @title Create an adaptor using indices
#'
#' @description
#' Creates an adaptor using a specified numeric vector of positive integer indices.
#' These indices are for subsetting operations.
#'
#' @param indices A numeric vector of positive integers representing the subsetting indices
#'
#' @return "adaptor" object
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' @seealso \code{\link{fixed_parameter_adaptor}}
indices_adaptor <- function(indices){
  if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("indices must be positive integers")
  adaptor(makeVectorSubsetByIndicesAdaptor(indices = indices))
}
