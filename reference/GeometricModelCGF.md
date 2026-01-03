# Create a Parametric Geometric CGF Object

Creates a CGF object for the Geometric distribution with success
probability \\prob(\theta)\\ given by a user-supplied function or
adaptor. The resulting object correspond to the random variable that
counts the number of failures before achieving the first success. It
supports i.i.d. and non-identical contexts with optional length
enforcement via `iidReps`.

## Usage

``` r
GeometricModelCGF(prob, iidReps = "any", ...)
```

## Arguments

- prob:

  A function (or adaptor) that accepts a single parameter vector `theta`
  and returns the success probability \\prob\\ (a scalar) or a vector of
  probabilities.

- iidReps:

  Either `"any"` (no forced dimension) or a positive integer specifying
  how many i.i.d. blocks are expected. Each block correspond to one copy
  of the geometric variables (or multiple if `prob` is a vector).

- ...:

  Additional arguments passed to the underlying CGF creation function

## Value

A `CGF` object
