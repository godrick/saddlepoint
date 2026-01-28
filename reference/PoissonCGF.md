# Poisson CGF object

Ready-to-use CGF for Poisson. Accepts scalar or vector `lambda`. If
`length(tvec)` is a multiple of `length(lambda)`, evaluation proceeds
(iidReps="any").

## Usage

``` r
PoissonCGF
```

## Format

An object of class `VectorizedFunctionsCGF` (inherits from `CGF`, `R6`)
of length 25.

## Examples

``` r
# Evaluate K at t = 0.1 for lambda = 2
# PoissonCGF$K(0.1, 2)
```
