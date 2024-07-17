##' @title BinomialCGF
##' @description
##' A CGF object for the binomial distribution that includes the CGF and its related derivatives.
##' When using this object, specify the `parameter_vector` as a vector with
##' the number of trials `n` and the probability of success `p` in that order.
##' 
##' @examples
##' \dontrun{
##'   # Compute the CGF and its first derivative at t = 0 for n = 10 trials and p = 0.5 
##'   cgf_value = BinomialCGF$K(tvec = 0, parameter_vector = c(10, 0.5))
##'   first_derivative = BinomialCGF$K1(tvec = 0, parameter_vector = c(10, 0.5))
##'   }
##' @export 
BinomialCGF <- NULL

##' @title PoissonCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for a binomial distribution
##' @examples
##' \dontrun{
##'   # Access the CGF for a Poisson distribution with lambda = 15
##'   cgf_value = PoissonCGF$K(tvec = 0.1, parameter_vector = 15)
##'   }
##' @export 
PoissonCGF <- NULL

##' @title ExponentialCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for the exponential distribution
##' @examples
##' \dontrun{
##'   mean_value = ExponentialCGF$K1(tvec = 0, parameter_vector = 15)
##'   }
##' @export 
ExponentialCGF <- NULL

##' @title GeometricCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for the geometric distribution
##' @examples
##' \dontrun{
##'   ## For X~Geometric(0.03), E(X) is  
##'   GeometricCGF$K1(tvec = 0, parameter_vector = 0.03) 
##'   }
##' @export 
GeometricCGF <- NULL

##' @title GammaCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for the gamma distribution
##' @examples
##' \dontrun{
##'   ## For X~Geometric(shape = 7.5, rate = 2), E(X) is  
##'   GammaCGF$K1(tvec = 0, parameter_vector = c(7.5, 2)) 
##'   }
##' @export 
GammaCGF <- NULL





.onLoad <- function(libname, pkgname) {
  utils::assignInMyNamespace('BinomialCGF', createCGF(make_BinomialCGF()))
  utils::assignInMyNamespace('PoissonCGF', createCGF(make_PoissonCGF()))
  utils::assignInMyNamespace('ExponentialCGF', createCGF(make_ExponentialCGF()))
  utils::assignInMyNamespace('GeometricCGF', createCGF(make_GeometricCGF()))
  utils::assignInMyNamespace('GammaCGF', createCGF(make_GammaCGF()))
}
