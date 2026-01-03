# Create a Parametric Binomial CGF Object

Binomial with parameters `n(theta)` and `prob(theta)`. If
`iidReps = "any"` (default), `length(tvec)` must be a multiple of
`length(n(theta))`. If `iidReps = m`, then `length(tvec)` must be
`m * length(n(theta))`. It supports i.i.d. replication (via `iidReps`)
and non-identical usage by causing `n(theta)` and `prob(theta)` to
return vectors.

## Usage

``` r
BinomialModelCGF(n, prob, iidReps = "any", ...)
```

## Arguments

- n:

  A function or `adaptor` mapping `theta` -\> scalar or vector of trial
  counts.

- prob:

  A function or `adaptor` mapping `theta` -\> scalar or vector of
  probabilities.

- iidReps:

  Either `"any"` or a positive integer. Default `"any"`.

- ...:

  Passed through to the CGF creation (e.g., to override operators).

## Value

A `CGF` object.

## Examples

``` r
# i.i.d. example: n=10, p=0.3
n_fn <- function(th) th[1]
p_fn <- function(th) th[2]
my_cgf <- BinomialModelCGF(n_fn, p_fn, iidReps = 1)
# If iidReps="any", length(tvec) determines the number of i.i.d. replicates
my_cgf$K1(0, c(10, 0.3))
#> [1] 3

# non-identical example: n=c(5,10), p=c(0.2,0.7)
# param ==> c(5,10, 0.2, 0.7)
n_adapt <- function(th) th[1:2] # OR n_adapt <- adaptor(indices = 1:2)
p_adapt <- function(th) th[3:4] # OR p_adapt <- adaptor(indices = 3:4)
my_cgf2 <- BinomialModelCGF(n_adapt, p_adapt, iidReps="any")
# default iidReps="any".
# If iidReps=m, then m-i.i.d. blocks are expected => length(tvec)
# must be a multiple of m
# e.g. tvec= c(0,0) => length=2 => we have 2 binomial distributions: (5,0.2) and (10,0.7)
my_cgf2$K1(c(0,0), c(5,10,0.2,0.7))
#> [1] 1 7
```
