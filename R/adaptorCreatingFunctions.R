#' @title Create an adaptor using fixed values
#'
#' @description
#' Creates an adaptor using a specified numeric vector of fixed parameter values.
#' This is especially useful for setting known values of parameters that will not change during estimation.
#' For instance, one might set the 'n' parameter to 1 for a Bernoulli distribution (which is a special case of the binomial distribution).
#'
#'
#' @param fixedParam A numeric vector indicating the fixed parameter values.
#'
#' @return An adaptor object, which is an S4 object of class "adaptor".
#'
#' @examples
#' \dontrun{
#'   # Define a fixed parameter for a Bernoulli trial
#'   fixedValue <- 1
#'   # Create an adaptor using the fixed parameter
#'   myAdaptor <- adaptorUsingFixedParam(fixedParam = fixedValue)
#'   class(myAdaptor) # should return "adaptor"
#'
#'   # How to use the adaptor in practice:
#'   # Here, we set the 'n' parameter of a Binomial distribution as a fixed value, and
#'   # the 'prob' parameter is set using an index pointing adaptor
#'   cgf = BinomialModelCGF(n = myAdaptor, prob = adaptor.using.indices(1))
#'   # Evaluate K using the created CGF object
#'   K(tvec = 0, parameter.vector = 0.5, baseCGF = cgf)
#' }
#'
#' @export
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}
adaptorUsingFixedParam <- function(fixedParam){
  if(!is.numeric(fixedParam)) stop("fixedParam is not a numeric vector")
  adaptor(makeSavedVectorAdaptor(fixed_parameter_values = fixedParam))
}



#' @title Create an adaptor using indices
#'
#' @description
#' Creates an adaptor using a specified numeric vector of positive integer indices.
#' These indices are intended for subsetting operations.
#' This function belongs to a suite of functions designed to generate adaptor objects.
#'
#' @param indices A numeric vector of positive integers representing the subsetting indices.
#'
#' @return An adaptor object, which is an S4 object of class "adaptor".
#'
#' @examples
#' \dontrun{
#'   # Create an adaptor using an index pointing to the third position in the parameter vector
#'   myAdaptor <- adaptorUsingIndices(indices = 3)
#'   class(myAdaptor) # should return "adaptor"
#'
#'   # Using the adaptors
#'   # Here, the 'n' parameter is set to a fixed value of 1 using (\code{\link{adaptorUsingFixedParam}})
#'   # and the 'prob' parameter will be fetched from the third position of the
#'   # `parameter.vector` due to the index specified in `myAdaptor`.
#'   cgf = BinomialModelCGF(n = adaptorUsingFixedParam(1), prob = myAdaptor)
#'
#'   # The parameter.vector has the value 0.5 at the third position, which will be used as the 'prob' parameter
#'   K(tvec = 0, parameter.vector = c(345, 9, 0.5), baseCGF = cgf)
#' }

#'
#' @export
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingFixedParam}}
adaptorUsingIndices <- function(indices){
  if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("indices must be positive integers")
  adaptor(makeVectorSubsetByIndicesAdaptor(indices = indices))
}


#' @title Create an adaptor Using R Functions
#'
#' @description
#' Creates an adaptor object using provided R functions. This is useful when specific transformations
#' of the parameter space are desired, or for implementing a derivative of the function with respect to parameters.
#'
#' @param h A function that computes h(theta). The function should accept a numeric vector as input and return a numeric vector.
#' @param grad_h A function that computes the gradient of h(theta) with respect to theta. It should accept a numeric vector as input and return a matrix with dimensions length(h(theta)) rows and length(theta) columns
#'
#' @return An adaptor object, which is an S4 object of class "adaptor".
#'
#' @examples
#' \dontrun{
#'   # Example functions h and grad.h
#'   h = function(theta) c(theta[1], 5, 3, theta[2], theta[3]) # theta is a numeric vector of length 3
#'   grad.h = function(theta) matrix(c(1,0,0,0,0, 0,0,0,1,0, 0,0,0,0,1), ncol = 3)
#'
#'   myAdaptor1 = adaptorUsingRFunctions(h, grad.h)
#'   class(myAdaptor1) # should return "adaptor"
#'
#'   # Another set of functions
#'   f = function(theta) c(theta^2, 2*theta, theta)
#'   grad.f = function(theta) {
#'      matrix(c(2*theta, 2, 1), ncol = 1)
#'   }
#'   myAdaptor2 = adaptorUsingRFunctions(f, grad.f)
#'   class(myAdaptor2) # should return "adaptor"
#' }
#'
#' @export
#' @seealso \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
adaptorUsingRFunctions <- function(h, grad_h){
  if(!is.function(h) || !is.function(grad_h)) stop("Both h and grad_h must be functions")
  adaptor(makeAdaptorUsingRfunctions(fn = h, grad_fn = grad_h))
}


