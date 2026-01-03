# Create a Parametric Exponential CGF Object

Constructs a CGF object for the Exponential distribution where the rate
parameter \\\lambda(\theta)\\ is derived from a user-supplied function
(or adaptor). This supports both i.i.d. and non-identical usage,
depending on whether \\\lambda(\theta)\\ returns one or multiple values,
and depending on the `iidReps` setting.

## Usage

``` r
ExponentialModelCGF(rate, iidReps = "any", ...)
```

## Arguments

- rate:

  A function (or adaptor) that accepts a parameter vector `theta` and
  returns the rate parameter \\\lambda\\ (a positive numeric value or
  vector).

- iidReps:

  Either `"any"` or a positive integer specifying the number of i.i.d.
  blocks are expected. Each block correspond to one copy of the
  Exponential variables (or multiple if \\\lambda(\theta)\\ is a
  vector).

- ...:

  Additional arguments passed to the underlying CGF creation function.

## Value

A `CGF` object.

## Examples

``` r
rate_func <- function(theta) theta[1]  # For example, theta -> 2 gives rate = 2
expo_model_cgf <- ExponentialModelCGF(rate = rate_func, iidReps = 1)
# Evaluate the first derivative at t = 0 for rate 2:
expo_model_cgf$K1(0, c(2))
#> [1] 0.5
```
