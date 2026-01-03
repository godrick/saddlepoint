# Create a parametric univariate Normal CGF

Builds a CGF for \\Y \sim N(\mu(\theta), \sigma(\theta))\\, where `mu`
and `sigma` are functions (or adaptors) of `theta`.

Both `mu(theta)` and `sigma(theta)` may return a scalar or a vector. If
they return vectors, they must have the same length L and are
interpreted in the same vectorised way as other univariate model CGFs in
this package.

## Usage

``` r
NormalModelCGF(mu, sigma, iidReps = "any", ...)

GaussianModelCGF(mu, sigma, iidReps = "any", ...)
```

## Arguments

- mu:

  Function/adaptor mapping theta -\> mean(s).

- sigma:

  Function/adaptor mapping theta -\> sd(s).

- iidReps:

  Either "any" or a positive integer (replication semantics).

- ...:

  Passed through to CGF creation.
