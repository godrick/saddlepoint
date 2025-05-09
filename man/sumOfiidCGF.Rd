% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/D-SumOfiidCGF.R
\name{sumOfiidCGF}
\alias{sumOfiidCGF}
\title{Create a CGF object for the sum of \code{n} i.i.d. random variables.}
\usage{
sumOfiidCGF(cgf, n, block_size = NULL, iidReps = NULL, ...)
}
\arguments{
\item{cgf}{A \code{CGF} object describing the base random variable \eqn{X}.}

\item{n}{Either a positive scalar or a function(\code{theta}) -> (positive scalar).}

\item{block_size}{Either \code{NULL} or a positive integer specifying the block size
for i.i.d. replication of the resulting CGF. Defaults to \code{NULL}.}

\item{iidReps}{Either \code{NULL} or a positive integer specifying how many i.i.d. blocks will be expected by the resulting CGF. Defaults to \code{NULL}.}

\item{...}{Additional arguments passed to the \code{CGF} creation function}
}
\value{
A new \code{CGF} object for \eqn{\;Y = \sum_{i=1}^{n(\theta)} X_i.}
}
\description{
Given a base CGF object \code{cgf} describing a random variable \eqn{X},
this function returns a new CGF object for the sum of \code{n} i.i.d. copies of \eqn{X},
i.e. \code{Y = X1 + X2 + ... + Xn}.

The argument \code{n} can be:
\itemize{
\item A positive numeric scalar (e.g., \code{n=5}).
\item A function or adaptor, \eqn{\theta \mapsto n(\theta)},
returning a positive numeric value for some parameter vector \code{theta}.
}
}
\examples{
\dontrun{
# Example: If 'cgf' describes a single Poisson(\lambda),
# then sumOfiidCGF(cgf, n=5) describes Poisson(5*lambda).

# Evaluate K, K1, etc. on the new object at t=0.1
# new_cgf <- sumOfiidCGF(cgf, n=5)
# new_cgf$K(0.1, param= c(...)) # => 5 times the old K(0.1)
}
}
