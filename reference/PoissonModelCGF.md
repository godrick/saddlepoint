# Create a Parametric Poisson CGF Object

Poisson with rate(s) `lambda(theta)`. If `iidReps = "any"` (default),
`length(tvec)` must be a multiple of `length(lambda(theta))`. Supports
i.i.d. and non-identical contexts with optional length enforcement via
`iidReps`. If `iidReps = m`, then `length(tvec)` must be
`m * length(lambda(theta))`.

## Usage

``` r
PoissonModelCGF(lambda, iidReps = "any", ...)
```

## Arguments

- lambda:

  A function or `adaptor` mapping `theta` -\> scalar or vector of rates.

- iidReps:

  Either `"any"` or a positive integer. Default `"any"`.

- ...:

  Passed through to the CGF creation (e.g., to override operators).

## Value

A `CGF` object.

## Examples

``` r
# ex 1: i.i.d. scenario for univariate Poisson(rate = 2)
ex_iid <- PoissonModelCGF(lambda = function(x) x[1])
# OR ex_iid <- PoissonModelCGF(lambda = function(x) x[1], iidReps = 3)
ex_iid$K1(rep(0,3), 2)
#> [1] 2 2 2

# ex 2: non-identical scenario with replication: lambda returns c(2,5), iidReps=3
ex_repeat <- PoissonModelCGF(lambda = adaptor(indices = 1:2), iidReps = 3)
# OR ex_repeat <- PoissonModelCGF(lambda = function(x) c(2,5) )
ex_repeat$K1(rep(0,6), c(2,5))
#> [1] 2 5 2 5 2 5
```
