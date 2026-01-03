# Compute the correction term for the saddlepoint approximation to the log-likelihood

This function calculates the correction term for the saddlepoint
approximation to the log-likelihood. By default, it computes that of the
**first-order** approximation, If `spa_method="zeroth"`, it instead
computes the correction term for the **zeroth-order** approximation.

## Usage

``` r
compute.saddlepointLL.correction(
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

  Numeric vector of observed data.

- cgf:

  A `CGF` object.

- tvec.hat:

  (Optional) Numeric vector. If supplied, `tvec` is taken directly as
  this vector (the saddlepoint \\\hat{t}\\). Otherwise we compute it
  (analytic if available, else numeric solve).

- gradient:

  Logical. If `TRUE`, return gradient wrt `parameter_vector`.

- hessian:

  Logical. If `TRUE`, return Hessian wrt `parameter_vector`.

- spa_method:

  Character string. One of `"standard"` (first-order) or `"zeroth"`.

- tvec_source:

  One of "auto","user","analytic","newton","solver_atomic". See
  `compute.spa.negll`.

- solver_fun:

  Numeric solver used when `tvec_source="solver_atomic"` or in non-AD
  path. Default `saddlepoint.solve`.

- newton_t_init:

  Optional start for t in the Newton solver path.

- ...:

  Passed to `solver_fun` when used.

## Value

A named list with elements:

- val:

  Numeric scalar of the correction term.

- gradient:

  The gradient at `theta`, if `gradient=TRUE`.

- hessian:

  The Hessian at `theta`, if `hessian=TRUE`.

## See also

[`find.saddlepoint.MLE`](https://godrick.github.io/saddlepoint/reference/find.saddlepoint.MLE.md),
[`compute.spa.negll`](https://godrick.github.io/saddlepoint/reference/compute.spa.negll.md)
