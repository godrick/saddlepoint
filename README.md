# saddlepoint

## Overview

The `saddlepoint` package (v0.9.0) provides tools for working with Cumulant Generating Functions (CGF) in statistical modeling, focusing on saddlepoint approximation methods. It's designed to automate computation of CGFs, even composite random variables. For parameter estimations using saddlepoint likelihood, automatic differentiation tool is used for accurate gradient-based optimizations.

**Notice:** This version is under active development. Significant changes to features and functionality are expected.

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
K(tvec = 0, parameter.vector = c(10, 0.5), baseCGF = GammaCGF)

# First derivative (mean)
K1(tvec = 0, parameter.vector = c(10, 0.5), baseCGF = GammaCGF)

# Second derivative (variance)
K2(tvec = 0, parameter.vector = c(10, 0.5), baseCGF = GammaCGF)

# Sample data from Gamma distribution
x = rgamma(50, shape = 10, rate = 0.5)

# MLE using saddlepoint likelihood
find.saddlepoint.MLE(observed.data = x, model.cgf = GammaCGF, starting.theta = c(1,1))$MLEs.theta


