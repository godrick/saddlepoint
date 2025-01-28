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


##' @title NegBinCGF
##' @description
##' A CGF object for the negative binomial distribution. 
##' When using this object, specify the `parameter_vector` as a vector with the number of successes `r` and the probability of success `p` in that order.
##' 
##' @export 
NegBinCGF <- NULL


##' @title PoissonCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for a Poisson distribution
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

##' @title MultinomialCGF
##' @description
##' An CGF object containing the CGF and its related derivatives for the multinomial distribution
##' @examples
##' \dontrun{
##'   ## For X~Multinomial(n = 10, prob = c(0.1, 0.2, 0.3, 0.4)), E(X) is  
##'   MultinomialCGF$K1(tvec = rep(0,4), parameter_vector = c(10, 0.1, 0.2, 0.3, 0.4)) 
##'   }
##' @export 
MultinomialCGF <- NULL






.onLoad <- function(libname, pkgname) {
  utils::assignInMyNamespace('BinomialCGF', createCGF(make_BinomialCGF()))
  utils::assignInMyNamespace('PoissonCGF', createCGF(make_PoissonCGF()))
  utils::assignInMyNamespace('ExponentialCGF', createCGF(make_ExponentialCGF()))
  utils::assignInMyNamespace('GeometricCGF', createCGF(make_GeometricCGF()))
  utils::assignInMyNamespace('GammaCGF', createCGF(make_GammaCGF()))
  utils::assignInMyNamespace('MultinomialCGF', createCGF(make_MultinomialCGF()))
  utils::assignInMyNamespace('NegBinCGF', createCGF(make_NegBinCGF()))
  
}
