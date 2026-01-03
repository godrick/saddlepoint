# Univariate Normal CGF object

Ready-to-use CGF for a univariate Normal distribution with parameters
\\(\mu,\sigma)\\.

## Usage

``` r
NormalCGF

GaussianCGF
```

## Format

An object of class `VectorizedFunctionsCGF` (inherits from `CGF`, `R6`)
of length 24.

An object of class `VectorizedFunctionsCGF` (inherits from `CGF`, `R6`)
of length 24.

## Details

The parameter vector is interpreted as `c(mu, sigma)` (length 2), or in
vectorised form as `c(mu[1:L], sigma[1:L])`.

With `iidReps="any"` (the default), `length(tvec)` must be a multiple of
the number of parameter rows (L). In the scalar case (L=1), any length
`tvec` is treated as i.i.d. replicates.
