# saddlepoint

## Overview

The `saddlepoint` package (v0.9.0) provides tools for working with Cumulant Generating Functions (CGF) in statistical modeling, focusing on saddlepoint approximation methods. It's designed to automate computation of CGFs, even composite random variables. For parameter estimations using saddlepoint likelihood, automatic differentiation tool is used for accurate gradient-based optimizations.



## Installation

```R
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("godrick/saddlepoint")
```


### Basic Usage Example with Gamma Distribution

This example demonstrates using the `saddlepoint` package to calculate the CGF, its first and second derivatives, and obtaining MLEs using the saddlepoint likelihood for a Gamma distribution.

#### CGF, derivatives and saddlepoint MLE
```R
# CGF at t = 0, shape = 10, rate = 0.5
GammaCGF$K(tvec = 0, parameter_vector = c(10, 0.5))

# First derivative (mean)
GammaCGF$K1(tvec = 0, parameter_vector = c(10, 0.5))

# Second derivative (variance)
GammaCGF$K2(tvec = 0, parameter_vector = c(10, 0.5))

# Sample data from Gamma distribution
set.seed(1); x = rgamma(50, shape = 10, rate = 0.5)

# MLE using saddlepoint likelihood
find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF, starting.theta = c(1,1))$MLEs.theta


