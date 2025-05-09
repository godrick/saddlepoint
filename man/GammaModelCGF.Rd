% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C-GammaCGF.R
\name{GammaModelCGF}
\alias{GammaModelCGF}
\title{Create a Parametric Gamma CGF Object}
\usage{
GammaModelCGF(shape, rate, iidReps = "any", ...)
}
\arguments{
\item{shape}{A function (or \code{adaptor})  that accepts a single parameter vector \code{theta} and returns the shape parameter.}

\item{rate}{A function (or \code{adaptor})  that accepts a single parameter vector \code{theta} and returns the rate parameter.}

\item{iidReps}{Either \code{"any"} or a positive integer specifying how many
i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.}

\item{...}{Additional arguments passed to the underlying CGF creation function.}
}
\value{
A \code{CGF} object.
}
\description{
Creates a CGF object for the Gamma distribution with shape \eqn{\alpha(\theta)} and
rate \eqn{\beta(\theta)} defined by user-provided parameter functions.
This function supports both i.i.d. and non-identical usage.
}
