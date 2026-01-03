# Create a Parametric MultinomialCGF Object

Constructs a CGF object for the multinomial distribution based on
user-specified functions or adaptors that map a parameter vector `theta`
to the multinomial's total count \\n\\ and probability vector \\p\\.
Optionally, you can replicate this CGF for `iidReps` i.i.d. multinomial
observations, thus imposing a length restriction on `tvec`.

## Usage

``` r
MultinomialModelCGF(n, prob_vec, iidReps = "any", ...)
```

## Arguments

- n:

  An adaptor or function defining the total count \\n(\theta)\\. It must
  accept a parameter vector \\\theta\\ and return a single numeric
  value.

- prob_vec:

  An adaptor or function defining the probability (or odds) vector
  \\p(\theta)\\. It must accept a parameter vector \\\theta\\ and return
  a numeric vector.

- iidReps:

  Either `"any"` or a positive integer specifying how many

- ...:

  Additional named arguments passed to
  [`createCGF`](https://godrick.github.io/saddlepoint/reference/createCGF.md),
  such as method overrides or operator definitions.

## Value

A `CGF` object (an R6 class) specialized for the multinomial
distribution with user-defined mappings. The returned object supports
the usual CGF methods (`K`, `K1`, `K2`, `K3operator`, etc.).

## Details

**User-Specified Parameter Mapping**:

- `n`: An adaptor or a function of the form
  `function(theta) -> numeric`. Must return a single numeric value,
  interpreted as the total count \\n\\ for the multinomial distribution.
  The returned value should be non-negative.

- `prob_vec`: An adaptor or a function of the form
  `function(theta) -> numeric vector`. Must return a vector of
  non-negative entries, interpreted as probabilities or odds for the
  multinomial distribution. If \\prob_vec\\ sums to 1, it is treated as
  a probability vector; otherwise, it is normalized into probabilities.

**I.I.D. Replicates**: By setting `iidReps` to a positive integer \\m\\,
you declare that the input vector \\tvec\\ will be split into \\m\\
blocks of equal size, each corresponding to one i.i.d. multinomial
sample. If `iidReps` is `"any"`, the number of blocks is inferred from
the input. It requires \\\mathrm{length}(tvec)\\ to be a multiple of \\d
= \mathrm{length}(prob\\vec(\theta))\\ (the category dimension), and
then sets \\m = \mathrm{length}(tvec)/d\\.

## Examples

``` r
# Suppose we want n = 10, and probabilities = c(0.2, 0.3, 0.5).
n_adaptor <- adaptor(fixed_param = 10)
p_adaptor <- adaptor(fixed_param = c(0.2, 0.3, 0.5))

my_cgf <- MultinomialModelCGF(n = n_adaptor, prob_vec = p_adaptor, iidReps = 2)

# The resulting CGF object expects 2 i.i.d. multinomial blocks,
# each of dimension 3. So tvec must be of length 6.
tvec <- rep(0, 6)
param <- c(10, 0.2, 0.3, 0.5)  # this will not be used since we have fixed_param

# We can now compute, e.g., the gradient:
my_cgf$K1(tvec, param)
#> [1] 2 3 5 2 3 5

```
