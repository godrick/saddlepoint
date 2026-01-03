# Geometric CGF Object

A ready-to-use CGF object for a single-parameter Geometric distribution.
This corresponds to the count of failures before the first success. By
default, this object is vectorized for i.i.d. replicates of probability
`prob`.

## Usage

``` r
GeometricCGF
```

## Format

An object of class `CGF` (R6), with usual methods:
`K, K1, K2, K3operator, K4operator`, etc.

## See also

[`GeometricModelCGF`](https://godrick.github.io/saddlepoint/reference/GeometricModelCGF.md)

## Examples

``` r
# expected value of X~Geometric(prob = 0.3) via CGF
GeometricCGF$K1(0, 0.3)
#> [1] 2.333333
```
