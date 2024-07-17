
validate_and_transform_adaptor <- function(obj) {
  if (is.function(obj)) {
    if (length(formals(obj)) != 1) {
      stop("The function must have exactly one argument.")
    }
    return(makeAdaptorUsingRfunctions(obj))
  }
  if (!inherits(obj, "adaptor")) stop("Provided argument must be either a function or an object of type 'adaptor'.")
  obj
}


#' @title BinomialCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the binomial distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic parameter adjustments in three ways: 
#' 1. R functions for dynamically computing parameters.
#' 2. adaptor objects using indices to specify parameters (\code{\link{indices_adaptor}}).
#' 3. adaptor objects that set fixed values (\code{\link{fixed_parameter_adaptor}}).
#' For static parameter values, it is recommended to use `fixed_parameter_adaptor`.
#'
#'
#'
#' @param n adaptor or function representing the number of trials. Automatically converted
#'          to an adaptor if provided as a function.
#' @param prob adaptor or function representing the probability parameter. Automatically converted
#'          to an adaptor if provided as a function.
#' @return 'CGF' object.
#'
#' @export
#' @seealso \code{\link{BinomialCGF}}
#'
#' @examples
#' \dontrun{
#' ## Example using a fixed parameter for `n` and an index adaptor for `prob`
#' # `fixed_parameter_adaptor` sets `n` to a constant value of 10
#' n_adaptor <- fixed_parameter_adaptor(10)  
#' 
#' # `indices_adaptor` assumes `prob` is part of a parameter vector that will be accessed using the provided index
#' # For example, if the parameter vector at runtime is c(0.5), then `indices_adaptor(1)` will use the first element.
#' prob_adaptor <- indices_adaptor(1)  
#' 
#' # Create the CGF object with these adaptors
#' cgf <- BinomialModelCGF(n = n_adaptor, prob = prob_adaptor)
#' 
#' # Compute the first derivative of the CGF at tvec = 0, assuming the actual value of `prob` is 0.2
#' cgf$K1(tvec = 0, parameter_vector = 0.2)
#' 
#' ## Example using a fixed parameter for `n` and a function for `prob`
#' # Define a function for `prob` that extracts the third element from a given vector
#' # This assumes that the parameter vector passed to the function will always have at least three elements
#' prob_function <- function(x) x[3] # returns the third element, which is expected to be the probability
#' 
#' # Create the CGF object
#' cgf1 <- BinomialModelCGF(n = n_adaptor, prob = prob_function)
#' 
#' # Example usage, assuming the parameter vector is supplied at runtime
#' # For instance, the parameter vector might be supplied during a calculation call
#' example_params <- c(0.1, 0.7, 0.2)  # Example vector, where 0.2 is the probability of success
#' cgf1$K1(tvec = 0, parameter_vector = example_params)
#' 
#'}
BinomialModelCGF <- function(n, prob){
  n = validate_and_transform_adaptor(n)
  prob = validate_and_transform_adaptor(prob)
  createCGF(make_BinomialModelCGF(n_adaptor = n, prob_adaptor = prob))
}







#' @title PoissonCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the Poisson distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic adjustments of the rate parameter (\code{lambda}) using similar mechanisms as described for the binomial distribution.
#' For detailed explanations on using adaptors for dynamic parameter adjustments, see \code{\link{BinomialModelCGF}}.
#'
#' @param lambda Adaptor or function representing the rate parameter of the Poisson distribution. 
#'               Automatically converted to an adaptor if provided as a function.
#' @return 'CGF' object configured for the Poisson distribution.
#'
#' @export
#' @seealso \code{\link{PoissonCGF}} \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#' # Example: lambda as a function
#' f <- function(x) { 2 + 0.1 * x }
#' cgf <- PoissonModelCGF(lambda = f)
#' cgf$K1(tvec = 0, parameter_vector = 3) == PoissonCGF$K1(tvec = 0, parameter_vector = 2 + 0.1 * 3)
#'}
PoissonModelCGF <- function(lambda) {
  lambda = validate_and_transform_adaptor(lambda)
  createCGF(make_PoissonModelCGF(lambda_adaptor = lambda))
}










#' @title ExponentialCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the exponential distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic adjustments of the rate parameter (\code{lambda}) using similar mechanisms as described for the binomial distribution.
#' For detailed explanations on using adaptors for dynamic parameter adjustments, see \code{\link{BinomialModelCGF}}.
#'
#' @param lambda Adaptor or function representing the rate parameter of the exponential distribution. 
#'               Automatically converted to an adaptor if provided as a function.
#' @return 'CGF' object configured for the exponential distribution.
#'
#' @export
#' @seealso \code{\link{ExponentialCGF}} \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#' # Example: lambda as a function
#' f <- function(x) { 0.1 * x }
#' cgf <- ExponentialModelCGF(lambda = f)
#' cgf$K1(tvec = 0, parameter_vector = 5) == ExponentialCGF$K1(tvec = 0, parameter_vector = 0.1 * 5)
#'}
ExponentialModelCGF <- function(lambda) {
  lambda = validate_and_transform_adaptor(lambda)
  createCGF(make_ExponentialModelCGF(lambda_adaptor = lambda))
}








#' @title GeometricCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the geometric distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic adjustments of the probability of success (\code{prob}) using similar mechanisms as described for the binomial distribution.
#' For detailed explanations on using adaptors for dynamic parameter adjustments, see \code{\link{BinomialModelCGF}}.
#'
#' @param prob Adaptor or function representing the probability of success parameter of the geometric distribution. 
#'               Automatically converted to an adaptor if provided as a function.
#' @return 'CGF' object configured for the geometric distribution.
#'
#' @export
#' @seealso \code{\link{GeometricCGF}} \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#' # Example: probability of success as a function
#' f <- function(x) { log(x) }
#' cgf <- GeometricModelCGF(prob = f)
#' cgf$K1(tvec = 0, parameter_vector = 1.05) == GeometricCGF$K1(tvec = 0, parameter_vector = log(1.05))
#'}
GeometricModelCGF <- function(prob) {
  prob = validate_and_transform_adaptor(prob)
  createCGF(make_GeometricModelCGF(p_adaptor = prob))
}




#' @title GammaCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the gamma distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic adjustments of the (\code{shape}) and (\code{rate}) parameters using similar mechanisms as described for the binomial distribution.
#' For detailed explanations on using adaptors for dynamic parameter adjustments, see \code{\link{BinomialModelCGF}}.
#'
#' @param shape Adaptor or function representing the shape parameter of the gamma distribution. 
#' @param rate Adaptor or function representing the rate parameter of the gamma distribution.
#'        
#' @return 'CGF' object configured for the gamma distribution.
#'
#' @export
#' @seealso \code{\link{GammaCGF}} \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#' ## Example using a fixed parameter for `rate` and an index adaptor for `shape`
#' # `rate` parameter fixed at 1
#' # `shape` will be part of a parameter vector that will be accessed using the the first index
#' cgf <- GammaModelCGF(shape = indices_adaptor(1), rate = fixed_parameter_adaptor(1))
#' cgf$K1(tvec = 0, parameter_vector = 0.2) == GammaCGF$K1(tvec = 0, parameter_vector = c(0.2, 1))
#'}
GammaModelCGF <- function(shape, rate) {
  shape_ = validate_and_transform_adaptor(shape)
  rate_ = validate_and_transform_adaptor(rate)
  createCGF(make_GammaModelCGF(shape_adaptor = shape_, rate_adaptor = rate_))
}


