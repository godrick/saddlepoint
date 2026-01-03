# Binomial CGF Object

Ready-to-use CGF for Binomial(n, p). Accepts packed params
`c(n_vec, p_vec)`. If `length(tvec)` is a multiple of `length(n_vec)`,
evaluation proceeds (iidReps = "any").

## Usage

``` r
BinomialCGF
```

## Format

An object of class `CGF` (R6), with standard methods: `K`, `K1`, `K2`,
`K3operator`, etc.

## Details

By default, `BinomialCGF` supports vectorized evaluation for i.i.d.
replicates.

## Examples

``` r
# Expected value of X ~ Binomial(10, 0.3)
BinomialCGF$K1(0, c(10, 0.3))
#> [1] 3
```
