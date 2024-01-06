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

