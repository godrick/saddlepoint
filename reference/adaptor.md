# Create an `adaptor` object

This function constructs an `adaptor` object from exactly one of:

1.  `indices`: A numeric vector of positive integer indices, used to
    subset `theta`.

2.  `fixed_param`: A numeric value (or vector) to be returned regardless
    of `theta`.

3.  `r_func`: A custom R function of the form
    `function(theta) -> numeric`.

## Usage

``` r
adaptor(indices = NULL, fixed_param = NULL, r_func = NULL)
```

## Arguments

- indices:

  A numeric vector of positive integers.

- fixed_param:

  A numeric vector or scalar to return.

- r_func:

  A function(`theta`) -\> numeric.

## Value

An object of class `adaptor`.

## Details

Exactly **one** of `indices`, `fixed_param`, or `r_func` must be
non-`NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example 1: Fixed and index-based parameter adaptors
  # Create a CGF object for a binomial distribution with a fixed parameter `n = 10`.
  # This configuration allows passing only the `prob` parameter in subsequent uses.
  binom_cgf <- BinomialModelCGF(n = adaptor(fixed_param = 10), prob = adaptor(indices = 1))
  # The object binom_cgf will expect only `prob` to be passed,
  # since `n` is already set.
  binom_cgf$K1(tvec = 0, parameter_vector = 0.3)

  # Example 2: Function adaptor
  # Dynamically compute probability parameter based on the input vector.
  cgf <- BinomialModelCGF(n = adaptor(indices = 2), prob = function(y) y[1])
  cgf$K1(tvec = 0, parameter_vector = c(0.3, 10))
} # }
```
