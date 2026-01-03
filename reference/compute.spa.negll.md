# Compute the saddlepoint negative log-likelihood

Computes either the "standard" or "zeroth-order" saddlepoint
approximation to the negative log–likelihood, with options for how the
saddlepoint t-vector is obtained.

## Usage

``` r
compute.spa.negll(
  parameter_vector,
  observed.data,
  cgf,
  tvec.hat = NULL,
  gradient = FALSE,
  hessian = FALSE,
  spa_method = "standard",
  tvec_source = c("auto", "user", "analytic", "newton", "solver_atomic"),
  solver_fun = saddlepoint.solve,
  newton_t_init = NULL,
  ...
)
```

## Arguments

- parameter_vector:

  Numeric vector of parameters.

- observed.data:

  Numeric vector of observations.

- cgf:

  A `CGF` object.

- tvec.hat:

  Optional numeric t-vector (one-time tape).

- gradient, hessian:

  Logical flags for derivatives.

- spa_method:

  "standard" or "zeroth".

- tvec_source:

  One of "auto","user","analytic","newton","solver_atomic".

- solver_fun:

  Numeric solver used when tvec_source="solver_atomic". Default
  `saddlepoint.solve`.

- newton_t_init:

  Optional start for t in the Newton solver path.

- ...:

  Passed to solver_fun (when used).

## Value

list(vals=..., gradient=..., hessian=...)
