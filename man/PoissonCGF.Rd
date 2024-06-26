% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributionCGFs.R
\docType{data}
\name{PoissonCGF}
\alias{PoissonCGF}
\title{PoissonCGF}
\format{
An object of class \code{"CGF"}.
}
\usage{
PoissonCGF
}
\description{
A CGF object based on the Poisson distribution.
}
\examples{
\dontrun{
  # Using the PoissonCGF in the K function
  K(tvec = 0, parameter.vector = 15, baseCGF = PoissonCGF)
}

}
\keyword{datasets}
