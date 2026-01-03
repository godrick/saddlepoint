# Negative Binomial CGF Object

A ready-to-use CGF object for the Negative Binomial distribution with
number of successes \\r\\ and success probability \\p\\. The
`parameter_vector` is \\c(r, p)\\, and the actual CGF is \$\$K(t) =
r\\\[\log p - \log(1 - (1-p)\\e^t)\], \quad t \< -\log(1-p).\$\$

## Usage

``` r
NegBinCGF
```

## Format

An object of class `CGF` (R6) with methods `K`, `K1`, `K2`,
`K3operator`, `K4operator`, etc.

## Examples

``` r
NegBinCGF$K1(0, c(10, 0.25))  # E[X] = r(1-p)/p = 10 * 0.75 / 0.25 = 30
#> [1] 30
```
