# Create a Parametric Gamma CGF Object

Creates a CGF object for the Gamma distribution with shape
\\\alpha(\theta)\\ and rate \\\beta(\theta)\\ defined by user-provided
parameter functions. This function supports both i.i.d. and
non-identical usage.

## Usage

``` r
GammaModelCGF(shape, rate, iidReps = "any", ...)
```

## Arguments

- shape:

  A function (or `adaptor`) that accepts a single parameter vector
  `theta` and returns the shape parameter.

- rate:

  A function (or `adaptor`) that accepts a single parameter vector
  `theta` and returns the rate parameter.

- iidReps:

  Either `"any"` or a positive integer specifying how many i.i.d. blocks
  are expected. Defaults to `"any"`, meaning no restriction on the
  length of `tvec`.

- ...:

  Additional arguments passed to the underlying CGF creation function.

## Value

A `CGF` object.
