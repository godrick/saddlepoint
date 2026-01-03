# Numerical saddlepoint equation solver

Solves the saddlepoint equation \\K1(t; \theta) = y\\ for a given CGF
object. It uses numeric optimization (via `nloptr`) to find the vector
\\t\\ that satisfies \\K1(t; \theta) = y\\.

If the CGF object offers an analytical solution (see
`cgf_object$analytic_tvec_hat` if available), you may prefer that
approach instead.

## Usage

``` r
saddlepoint.solve(
  theta,
  y,
  cgf,
  starting.tvec = rep(0, times = length(y)),
  lb = rep(-Inf, times = length(y)),
  ub = rep(Inf, times = length(y)),
  sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1000, xtol_rel = 1e-12, print_level = 0),
  warn_residual = TRUE,
  tol = 1e-04
)
```

## Arguments

- theta:

  Numeric vector of model parameters for which the CGF is defined.

- y:

  Numeric vector of observations.

- cgf:

  A CGF object (class `"CGF"`).

- starting.tvec:

  Numeric start values for \\t\\, defaults to `0` for each element.

- lb, ub:

  Numeric vectors for lower and upper bounds of \\t\\, each having the
  same length as `starting.tvec`. Defaults are `-Inf` and `Inf`.

- sadd.eqn.opts:

  List of options for the `nloptr` optimizer. Defaults to
  `list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)`,
  with algorithm fixed at `"NLOPT_LD_SLSQP"`.

- warn_residual:

  Logical. If `TRUE`, compute the residual
  \\\\\mathit{K1}(t^\*;\theta)-y\\\_\infty\\ and warn if it exceeds
  `tol`.

- tol:

  Numeric tolerance for residual checks (default `1e-4`).

## Value

A numeric vector \\t\\ solving \\K1(t, \theta) = y\\.

## Details

Minimizes the function \\K(t; \theta) - \sum(t_i y_i)\\ to enforce
\\K1(t, \theta) = y\\.

## References

- [NLOpt Documentation](https://nlopt.readthedocs.io/en/latest/)

## Examples

``` r
if (FALSE) { # \dontrun{
# TO DO: Add examples
} # }
```
