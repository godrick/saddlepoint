# Gamma CGF Object

A ready-to-use CGF object for the Gamma distribution with shape
\\\alpha\\ and rate \\\beta\\. The `parameter_vector` used when calling
methods such as `K(tvec, parameter_vector)` should be a numeric vector
\\c(\alpha, \beta)\\.

## Usage

``` r
GammaCGF
```

## Format

An object of class `CGF` (an R6 class), with the usual methods:
`K, K1, K2, K3operator, K4operator`, etc.

## Details

**CGF**: For a Gamma random variable \\X\\ with shape \\\alpha\\ and
rate \\\beta\\, the cumulant generating function is: \$\$K(t;\alpha,
\beta) = -\alpha \\\log \bigl(1 - t/\beta\bigr), \quad t \< \beta.\$\$

**Parameter Vector**: The `parameter_vector` is assumed to have the form
\\(\alpha, \beta)\\. You must ensure that \\t \< \beta\\ for valid
evaluations.

## Examples

``` r
# Evaluate K at t=0.5 for shape=2, rate=2 (thus t<2).
# param = c(2, 2)
GammaCGF$K(0.5, c(2,2))
#> [1] 0.5753641
```
