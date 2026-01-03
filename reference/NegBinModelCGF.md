# Create a Parametric Negative Binomial CGF Object

Creates a CGF object for the Negative Binomial distribution with number
of successes \\r(\theta)\\ and success probability parameter
\\p(\theta)\\ defined by user-provided parameter functions. This
function supports both i.i.d. and non-identical usage.

## Usage

``` r
NegBinModelCGF(r, p, iidReps = "any", ...)
```

## Arguments

- r:

  A function (or `adaptor`) that accepts a single parameter vector
  `theta` and returns the number of successes.

- p:

  A function (or `adaptor`) that accepts a single parameter vector
  `theta` and returns a scalar success probability or a vector of
  success probabilities.

- iidReps:

  Either `"any"` or a positive integer specifying how many i.i.d. blocks
  are expected. Defaults to `"any"`, meaning no restriction on the
  length of `tvec`.

- ...:

  Additional arguments passed to the underlying CGF creation function.

## Value

A `CGF` object.
