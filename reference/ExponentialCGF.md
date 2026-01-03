# Exponential CGF Object

A ready-to-use CGF object for the Exponential distribution. This object
is vectorized for i.i.d. replicates of the same rate.

## Usage

``` r
ExponentialCGF
```

## Format

An object of class `CGF` (R6), with usual methods `K`, `K1`, `K2`,
`K3operator`, `K4operator`, etc.

## Examples

``` r
ExponentialCGF$K1(0, 2)  # (expected value for Exp(rate=2) is 1/2)
#> [1] 0.5
```
