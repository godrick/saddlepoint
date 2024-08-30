

# This function allows direct use of R functions as adaptors
# We also use rtmb_Xtra environment when a user supplies an R function
validate_and_transform_adaptor <- function(obj) {
  if (!is.function(obj)) {
    if (inherits(obj, "adaptor")) {return(obj)}
    stop("Provided argument must be either a function or an object of type 'adaptor'.")
  }
  if (length(formals(obj)) != 1) stop("The function used must have exactly one argument.")
  r_func_adaptor <- function(x) {
    base::attach(rtmb_Xtra, length(search()), name = "rtmb_Xtra-AD-overloads", warn=FALSE)
    on.exit(base::detach("rtmb_Xtra-AD-overloads"))
    obj(x)
  }
  makeAdaptorUsingRfunctions(r_func_adaptor)
}


#' @title BinomialCGF with flexible parameterization
#'
#' @description
#' Constructs a CGF object for the binomial distribution, supporting flexible parameter specifications through adaptors.
#' 
#' @details
#' This function allows for dynamic parameter adjustments in three ways: 
#' 1. **Function-based adaptors**: Use R functions to dynamically compute parameters based on input vector.
#' 2. **Index-based adaptors**: Subset parameters from a vector using specified indices `adaptor(indices = ...)`.
#' 3. **Fixed Value adaptors**: Set parameters to constant values `adaptor(fixed_param = ...)`.
#' For fixed parameter values, it is recommended to use `adaptor(fixed_param = ...)` and not function-based adaptors.
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
#' ## Example using a fixed parameter for `n` and an index-based adaptor for `prob`
#' n_adaptor <- adaptor(fixed_param = 10) # Set `n` to a constant value of 10
#' 
#' # `indices` argument of the `adaptor()` function supposes that the target (here `prob`) is part of a parameter vector that will be accessed using the provided indices
#' # For example, if the parameter vector at runtime is c(0.2, 7, 4), then `adaptor(indices = 1)` will use the first element for the `prob`.
#' prob_adaptor <- adaptor(indices = 1) 
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
#' example_params <- c(0.1, 0.7, 0.2)  # Example vector, where 0.2 is the probability to be used
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
#' @param rate Adaptor or function representing the rate parameter of the exponential distribution. 
#'               Automatically converted to an adaptor if provided as a function.
#' @return 'CGF' object configured for the exponential distribution.
#'
#' @export
#' @seealso \code{\link{ExponentialCGF}} \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#' # Example: rate as a function
#' f <- function(x) { 0.1 * x }
#' cgf <- ExponentialModelCGF(rate = f)
#' cgf$K1(tvec = 0, parameter_vector = 5) == ExponentialCGF$K1(tvec = 0, parameter_vector = 0.1 * 5)
#'}makeAdaptorUsingRfunctions
ExponentialModelCGF <- function(rate) {
  rate_ = validate_and_transform_adaptor(rate)
  createCGF(make_ExponentialModelCGF(lambda_adaptor = rate_))
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
#' cgf <- GammaModelCGF(shape = adaptor(indices = 1), rate = adaptor(fixed_param = 1))
#' cgf$K1(tvec = 0, parameter_vector = 0.2) == GammaCGF$K1(tvec = 0, parameter_vector = c(0.2, 1))
#'}
GammaModelCGF <- function(shape, rate) {
  shape_ = validate_and_transform_adaptor(shape)
  rate_ = validate_and_transform_adaptor(rate)
  createCGF(make_GammaModelCGF(shape_adaptor = shape_, rate_adaptor = rate_))
}



#' @title MultinomialCGF with flexible parameterization
#' 
#' @description
#' Constructs a CGF object for the multinomial distribution, supporting flexible parameter specifications through adaptors.
#' This function allows for dynamic adjustments of (\code{n} and \code{prob_vec}) parameters using similar mechanisms as described for the binomial distribution.
#' For detailed explanations on using adaptors for dynamic parameter adjustments, see \code{\link{BinomialModelCGF}}.
#' 
#' @param n Adaptor or function representing the number of trials parameter of the multinomial distribution.
#' @param prob_vec Adaptor or function representing the probability vector parameter of the multinomial distribution.
#' 
#' @return 'CGF' object configured for the multinomial distribution.
#' 
#' @export
#' @seealso \code{\link{MultinomialCGF}} \code{\link{BinomialModelCGF}}
MultinomialModelCGF <- function(n, prob_vec) {
  n_ = validate_and_transform_adaptor(n)
  prob_vec_ = validate_and_transform_adaptor(prob_vec)
  createCGF(make_MultinomialModelCGF(n_, prob_vec_))
}





#' @title Adapt a CGF object
#' 
#' @description This function modifies a given CGF object by applying an adaptor or a function for dynamic parameterization.
#' It allows for custom behavior in parameter handling, allowing flexibility in parameter specification.
#' 
#'
#' @param cgf A 'CGF' object to be adapted
#' @param adaptor An 'adaptor' object or a function that modifies or provides parameters dynamically.
#'
#'
#' @return A CGF object.
#' @seealso \code{\link{adaptor}}, \code{\link{BinomialModelCGF}}
#'
#' @examples
#' \dontrun{
#'   ## Adapting the BinomialCGF using indices
#'   # Suppose the parameter_vector is of any length (> 2) and we want to use the first two elements as parameters of the underlying Binomial distribution
#'   adapted_binomial_cgf <- adaptCGF(cgf = BinomialCGF,
#'                                    adaptor = adaptor(indices = 1:2) )
#'   ## Adapting the BinomialCGF using a function
#'   adapted_binomial_cgf1 <- adaptCGF(cgf = BinomialCGF,
#'                                     adaptor = function(x) { x[1:2] } )
#'   res1 <- adapted_binomial_cgf$K1(tvec = 0, parameter_vector = c(10, 0.5, 90, 4, 3))
#'   res2 <- adapted_binomial_cgf1$K1(tvec = 0, parameter_vector = c(10, 0.5, 5))
#' }
#' @export
adaptCGF <- function(cgf, adaptor){
  if (!is(cgf, "CGF")) stop("'cgf' is not defined for ", class(cgf))
  param_adaptor = validate_and_transform_adaptor(adaptor)
  createCGF(adapt_CGF(cgf$get_ptr(), param_adaptor))
}



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
#' @param n An object of class 'adaptor' representing the number of trials in the Multinomial distribution.
#' @param prob_vec An object of class 'adaptor' representing the probabilities associated with the non-zero outcomes in the modified Multinomial distribution. The sum of these probabilities should be strictly less than one to validate the occurrence of the known zero outcome.
#'
#' @return A 'CGF' object.
#' @seealso \code{\link{MultinomialModelCGF}}
#'
#' @export
SubunitaryMultinomialModelCGF <- function(n, prob_vec){
  n_ = validate_and_transform_adaptor(n)
  prob_vec_ = validate_and_transform_adaptor(prob_vec)
  createCGF(make_SubunitaryMultinomialModelCGF(n_, prob_vec_))
}


