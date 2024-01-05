#' @title K (CGF)
#' @description Computes the Cumulant Genrating function, K for given t-vector, parameter vector, and a CGF object.
#'
#' @param tvec A numeric vector of t-values.
#' @param parameter.vector A numeric vector of parameters.
#' @param baseCGF A CGF type.
#'
#'
#' @return A numeric of length 1.
#'
#' @examples
#' \dontrun{
#'   K(tvec = c(0.003, 0.2), parameter.vector = c(10,0.5), baseCGF = BinomialCGF)
#' }
#' @export
K <- function(tvec, parameter.vector, baseCGF) {
  if (!is.numeric(tvec)) stop("'tvec' is not defined for ", class(tvec))
  if (!is.numeric(parameter.vector)) stop(" 'parameter.vector' is not defined for ", class(parameter.vector))
  if (!is(baseCGF, "CGF")) stop(" 'baseCGF' must be of class 'CGF'.")
  .Call('_saddlepoint_K', PACKAGE = 'saddlepoint', tvec, parameter.vector, baseCGF)
}
#' @title K1
#' @description Computes the tvec-gradient of K for given t-values, a parameter vector, and a CGF object.
#'
#' @param tvec A numeric vector of t-values.
#' @param parameter.vector A numeric vector of parameters.
#' @param baseCGF A CGF type.
#'
#'
#' @return A numeric vector whose length equals the length of 'tvec'.
#'
#' @examples
#' \dontrun{
#'   K1(tvec = c(0.003, 0.2), parameter.vector = c(10,0.5), baseCGF = BinomialCGF)
#' }
#' @export
K1 <- function(tvec, parameter.vector, baseCGF) {
  if (!is.numeric(tvec)) stop("'tvec' is not defined for ", class(tvec))
  if (!is.numeric(parameter.vector)) stop("'parameter.vector' is not defined for ", class(parameter.vector))
  if (!is(baseCGF, "CGF")) stop("'baseCGF' must be of class 'CGF'")
  .Call('_saddlepoint_K1', PACKAGE = 'saddlepoint', tvec, parameter.vector, baseCGF)
}

#' @title K2
#' @description Computes the Hessian matrix of the function K with respect to 'tvec', for given t-values, a parameter vector, and a CGF object.
#'
#' @param tvec A numeric vector of t-values.
#' @param parameter.vector A numeric vector of parameters.
#' @param baseCGF A CGF type.
#'
#'
#' @return A square matrix with dimensions equivalent to the length of 'tvec'.
#'
#' @examples
#' \dontrun{
#'   K2(tvec = c(0.003, 0.2), parameter.vector = c(10,0.5), baseCGF = BinomialCGF)
#' }
#' @export
K2 <- function(tvec, parameter.vector, baseCGF) {
  if (!is.numeric(tvec)) stop("'tvec' is not defined for ", class(tvec))
  if (!is.numeric(parameter.vector)) stop("'parameter.vector' is not defined for ", class(parameter.vector))
  if (!is(baseCGF, "CGF")) stop("'baseCGF' must be of class 'CGF'")
  .Call('_saddlepoint_K2', PACKAGE = 'saddlepoint', tvec, parameter.vector, baseCGF)
}

#' @title Compute tvec-related constraints for a CGF Object
#'
#' @description The function computes constraints for the argument 'tvec' under a given CGF object, expressed in a form of inequalities that should be less than 0.
#'
#' @param tvec A numeric vector representing 't' values.
#' @param parameter.vector A numeric vector of parameters.
#' @param baseCGF A 'CGF' object.
#'
#' @return A numeric vector of 'tvec' related constraints. An empty vector if no constraints on 'tvec'.
#'
#' @examples
#' \dontrun{
#' # Create a Geometric CGF object as a function of probability parameter 'p'.
#' # By definition, this CGF is valid for 't < -log(1 - p)'
#' # Fix 'p' at 0.5, so the domain of the CGF is constrained to 't < -log(0.5)'.
#'
#' # Compute constraints for the Geometric CGF object. The function 'ineq_constraint'
#' # will transform the domain restriction 't < -log(0.5)' into the equivalent form
#' # 't + log(0.5) < 0'. Hence, for each 't' value in 'tvec', we get an inequality
#' # 't + log(0.5)'.
#' ineq_constraint(tvec = c(0.1, 0.2), parameter.vector = 0.5, baseCGF = GeometricCGF)
#' }
#' @export
ineq_constraint <- function(tvec, parameter.vector, baseCGF) {
  if (!is.numeric(tvec)) stop("'tvec' is not defined for ", class(tvec))
  if (!is.numeric(parameter.vector)) stop("'parameter.vector' is not defined for ", class(parameter.vector))
  if (!is(baseCGF, "CGF")) stop("'baseCGF' must be of class 'CGF'")
  .Call('_saddlepoint_ineq_constraint', PACKAGE = 'saddlepoint', tvec, parameter.vector, baseCGF)
}
#' @title Compute the saddlepoint negative log-likelihood
#'
#' @description The function computes the saddlepoint negative log-likelihood (neg_ll) at given t-values, a parameter vector, and a CGF object. See ... for MLEs
#'
#' @param tvec A numeric vector of t-values at which the neg_ll is computed.
#' @param parameter.vector A numeric vector of parameters.
#' @param baseCGF A 'CGF' object which serves as the base for computing neg_ll.
#'
#' @return A numeric value, which is the computed saddlepoint negative log-likelihood.
#'
#' @examples
#' \dontrun{
#' # Create a CGF object
#' # Compute neg_ll for given t-values and parameters
#' }
#' @export
neg_ll <- function(tvec, parameter.vector, baseCGF) {
  if (!is.numeric(tvec)) stop("'tvec' is not defined for ", class(tvec))
  if (!is.numeric(parameter.vector)) stop("'parameter.vector' is not defined for ", class(parameter.vector))
  if (!is(baseCGF, "CGF")) stop("'baseCGF' must be of class 'CGF'")
  .Call('_saddlepoint_neg_ll', PACKAGE = 'saddlepoint', tvec, parameter.vector, baseCGF)
}


