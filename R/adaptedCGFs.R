#' @title BinomialModelCGF
#'
#' @description
#' Constructs a CGF for the Binomial distribution with flexible parameter specifications using adaptors. The adaptors offer multiple methods to specify parameters:
#' 1. Using indices from an external parameter vector (\code{\link{adaptorUsingIndices}}).
#' 2. Setting fixed values that won't change during model computations (\code{\link{adaptorUsingFixedParam}}).
#' 3. Employing R functions for dynamic computation of parameters (\code{\link{adaptorUsingRFunctions}}).
#'
#' This adaptability facilitates model constructions where specific parameters might be varied during computations while keeping others constant.
#'
#'
#' @param n An object of class 'adaptor', representing the parameter n in a Binomial distribution.
#' @param prob An object of class 'adaptor', representing the probability parameter in a Binomial distribution.
#'
#' @return An object of class 'CGF'.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Using adaptorUsingFixedParam for n and adaptorUsingIndices for prob
#'   cgf = BinomialModelCGF(n = adaptorUsingFixedParam(10), prob = adaptorUsingIndices(1))
#'   # Compute CGF with: K(tvec = ..., parameter.vector = 0.5, baseCGF = cgf)
#'}
#'
BinomialModelCGF <- function(n, prob){
  if (!is(n, "adaptor") || !is(prob, "adaptor")) stop(" 'n and prob' must be of type adaptor")
  CGF(make_BinomialModelCGF(n_adaptor = n, prob_adaptor = prob))
}

#' Poisson Model-Based CGF
#'
#' Generates a CGF based on the Poisson distribution.
#'
#' @param lambda An object of class 'adaptor', representing the rate parameter lambda in a Poisson distribution.
#'
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
PoissonModelCGF <- function(lambda) {
  if (!is(lambda, "adaptor")) stop(" 'lambda' is not defined for ", class(lambda))
  CGF(make_PoissonModelCGF(lambda_adaptor = lambda))
}

#' Exponential Model-Based CGF
#'
#' Generates a CGF based on the Exponential distribution.
#'
#' @param rate An object of class 'adaptor', representing the rate parameter in the Exponential distribution.
#'
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
ExponentialModelCGF <- function(rate){
  if (!is(rate, "adaptor")) stop(" 'rate' is not defined for ", class(rate))
  CGF(make_ExponentialModelCGF(lambda_adaptor = rate))
}

#' Geometric Model-Based CGF
#'
#' Generates a CGF based on the Geometric distribution.
#'
#' @param prob An object of class 'adaptor', representing the probability parameter in a Geometric distribution.
#'
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
GeometricModelCGF <- function(prob){
  if (!is(prob, "adaptor")) stop(" 'prob' is not defined for ", class(prob))
  CGF(make_GeometricModelCGF(p_adaptor = prob))
}

#' Geometric Non-Identical Model-Based CGF
#'
#' Generates a CGF based on a non-identical Geometric distribution.
#'
#' @param prob An object of class 'adaptor', representing the probability parameter in a Geometric distribution.
#'
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
GeometricNonIdenticalModelCGF <- function(prob){
  if (!is(prob, "adaptor")) stop(" 'prob' is not defined for ", class(prob))
  CGF(make_GeometricNonIdenticalModelCGF(p_adaptor = prob))
}

#' @title GammaModelCGF
#'
#' @description
#' Constructs a CGF for the Gamma distribution with parameter specifications using adaptors. The adaptors offer multiple methods to specify parameters:
#' 1. Using indices from an external parameter vector (\code{\link{adaptorUsingIndices}}).
#' 2. Setting fixed values that won't change during model computations (\code{\link{adaptorUsingFixedParam}}).
#' 3. Employing R functions for dynamic computation of parameters (\code{\link{adaptorUsingRFunctions}}).
#'
#'
#' @param shape An object of class 'adaptor', representing the parameter shape in a Gamma distribution.
#' @param rate An object of class 'adaptor', representing the rate parameter in a Gamma distribution.
#'
#' @return An object of class 'CGF'.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#'}
#'
GammaModelCGF <- function(shape, rate){
  if (!is(shape, "adaptor") || !is(rate, "adaptor")) stop(" 'shape and rate' must be of type adaptor")
  CGF(make_GammaModelCGF(shape_adaptor = shape, rate_adaptor = rate))
}

#' Multinomial Model-Based CGF
#'
#' Generates a CGF based on the Multinomial distribution.
#'
#' @param n An object of class 'adaptor'.
#' @param prob_vec An object of class 'adaptor', representing a vector of Multinomial probabilities.
#'
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
MultinomialModelCGF <- function(n, prob_vec){
  if(!is(n, "adaptor")) stop(" 'n' must be of type adaptor")
  if(!is(prob_vec, "adaptor")) stop(" 'prob_vec' must be of type adaptor")
  CGF(make_MultinomialModelCGF(n_adaptor = n, probVector_adaptor = prob_vec))
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
#' @return An object of class 'CGF'.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @export
SubunitaryMultinomialModelCGF <- function(n, prob_vec){
  if(!is(n, "adaptor")) stop(" 'n' must be of type adaptor")
  if(!is(prob_vec, "adaptor")) stop(" 'prob_vec' must be of type adaptor")
  CGF(make_SubunitaryMultinomialModelCGF(n_adaptor = n, probVector_adaptor = prob_vec))
}


#' @title adaptCGF
#' @description Adapts a CGF object using a specified adaptor.
#'
#' @param baseCGF A 'CGF' object to be adapted
#' @param adaptor An 'adaptor' object created by one of the adaptor creating functions.
#'
#' @details
#' The function create...
#'
#' @return A CGF object.
#' @seealso \code{\link{adaptorUsingRFunctions}}, \code{\link{adaptorUsingIndices}}, \code{\link{adaptorUsingFixedParam}}
#'
#' @examples
#' \dontrun{
#'   # Adapting the BinomialCGF with a fixed parameter adaptor
#'   adaptedCGF <- adapt.CGF(baseCGF = BinomialCGF,
#'                           adaptor = adaptorUsingFixedParam(c(10, 0.5)) )
#'
#'   # Creating a Binomial Model CGF with the same parameters
#'   modelCGF <- BinomialModelCGF(n = adaptorUsingFixedParam(10),
#'                                prob = adaptorUsingFixedParam(0.5))
#'
#'   # Testing if the adaptedCGF and modelCGF give the same result when used as baseCGF in function K
#'   res1 <- K(tvec = 0.003, parameter.vector = c(10,0.5), baseCGF = adaptedCGF)
#'   res2 <- K(tvec = 0.003, parameter.vector = c(10,0.5), baseCGF = modelCGF)
#'   print(all.equal(res1, res2)) # Should print TRUE
#' }
#' @export
adaptCGF <- function(baseCGF, adaptor){
  if (!is(baseCGF, "CGF")) stop("'baseCGF' is not defined for ", class(baseCGF))
  if (!is(adaptor, "adaptor")) stop(" 'adaptor' is not defined for ", class(adaptor))
  CGF(adapt_CGF(base_cgf = baseCGF, adaptor = adaptor))
}

