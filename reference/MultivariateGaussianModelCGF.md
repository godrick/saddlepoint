# Create a Multivariate Normal CGF Object

Creates a CGF for a \\d\\-dimensional Normal distribution with mean
`mu(theta)` and covariance `sigma(theta)`.

Replication:

- `iidReps` refers to i.i.d. replication of the \\d\\-vector
  observation.

- If `iidReps="any"` (default), `length(tvec)` must be a multiple of
  \\d\\.

- If `iidReps=m` (integer), `length(tvec)` must equal \\d\*m\\.

## Usage

``` r
MultivariateNormalModelCGF(mu, sigma, iidReps = "any", ...)

MultivariateGaussianModelCGF(mu, sigma, iidReps = "any", ...)
```

## Arguments

- mu:

  Function/adaptor mapping `theta` to a numeric vector of length \\d\\.

- sigma:

  Function/adaptor mapping `theta` to a \\d\times d\\ covariance matrix.

- iidReps:

  Either `"any"` or a positive integer.

- ...:

  Passed through to CGF creation.

## Value

A `CGF` object.
