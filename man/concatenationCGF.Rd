% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CGF-operations.R
\name{concatenationCGF}
\alias{concatenationCGF}
\title{Concatenation of CGF objects.}
\usage{
concatenationCGF(...)
}
\arguments{
\item{...}{A sequence of alternating CGF objects and integers. Each CGF object must be immediately followed by its corresponding vector length, which denotes the dimension of the underlying random variable.}
}
\value{
A CGF object.
}
\description{
Computes the CGF object for a concatenated vector of independent random variables, with each represented by its respective CGF object and a corresponding dimension. It facilitates the modeling of a concatenated random vector, where each segment of the vector is specified by its CGF object and a corresponding dimension.
}
