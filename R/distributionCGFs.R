#' @title BinomialCGF
#' @description
#' A CGF object based on the Binomial distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#'   # Using the BinomialCGF in the K function
#'   K(tvec = 0, parameter.vector = c(10, 0.5), baseCGF = BinomialCGF)
#' }
#'
#' @export
BinomialCGF <- NULL

#' @title GammaCGF
#' @description
#' A CGF object based on the Gamma distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
GammaCGF <- NULL

#' @title PoissonCGF
#' @description
#' A CGF object based on the Poisson distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#'   # Using the PoissonCGF in the K function
#'   K(tvec = 0, parameter.vector = 15, baseCGF = PoissonCGF)
#' }
#'
#' @export
PoissonCGF <- NULL

#' @title ExponentialCGF
#' @description
#' A CGF object based on the Exponential distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#'   # Using the ExponentialCGF in the K function
#'   K(tvec = 0, parameter.vector = 10, baseCGF = ExponentialCGF)
#' }
#'
#' @export
ExponentialCGF <- NULL

#' @title GeometricCGF
#' @description
#' A CGF object based on the Geometric distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#'   # Using the GeometricCGF in the K function
#'   K(tvec = 0, parameter.vector = 0.03, baseCGF = GeometricCGF) # Note that tvec must be less than -log(1-parameter.vector)
#' }
#'
#' @export
GeometricCGF <- NULL

#' @title MultinomialCGF
#' @description
#' A CGF object based on the Multinomial distribution.
#'
#' @format
#' An object of class \code{"CGF"}.
#'
#' @examples
#' \dontrun{
#'   # Using the MultinomialCGF in the K function
#'   K(tvec = rep(0,3), parameter.vector = c(10, 0.1, 0.2, 0.7), baseCGF = MultinomialCGF)
#' }
#'
#' @export
MultinomialCGF <- NULL


.onLoad <- function(libname, pkgname) {
  # Create CGFs within the .onLoad function
  utils::assignInMyNamespace('BinomialCGF', CGF(make_BinomialCGF()))
  utils::assignInMyNamespace('GammaCGF', CGF(make_GammaCGF()))
  utils::assignInMyNamespace('PoissonCGF', CGF(make_PoissonCGF()))
  utils::assignInMyNamespace('GeometricCGF', CGF(make_GeometricCGF()))
  utils::assignInMyNamespace('MultinomialCGF', CGF(make_MultinomialCGF()))
  utils::assignInMyNamespace('ExponentialCGF', CGF(make_ExponentialCGF()))
}
