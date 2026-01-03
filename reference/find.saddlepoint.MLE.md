# Find maximum likelihood estimates using the saddlepoint likelihood.

This function finds maximum likelihood estimates (MLEs) under a
saddlepoint approximation, using constrained or two-step optimization.

By default (`method="constrained"`) it performs the joint constrained
optimisation over \\(t,\theta)\\: \$\$ \hat t, \hat\theta =
\arg\min\_{(t,\theta)} \\ \mathrm{NLL}\_{\mathrm{SPA}}(t,\theta) \quad
\text{s.t.}\quad K_1(t;\theta) = y, \$\$ using `NLOPT_LD_SLSQP`.

Optionally, you may set `method="two_step"` to optimise over \\\theta\\
only, where \\\hat t(\theta)\\ is computed internally for each
\\\theta\\ using the RTMB/Newton pathway already implemented in
`compute.spa.negll`. This typically reduces the outer problem dimension
and MAY be much faster.

## Usage

``` r
find.saddlepoint.MLE(
  observed.data,
  cgf,
  starting.theta,
  lb.theta = rep(-Inf, times = length(starting.theta)),
  ub.theta = rep(Inf, times = length(starting.theta)),
  starting.tvec = rep(0, times = length(as.numeric(observed.data))),
  lb.tvec = rep(-Inf, times = length(as.numeric(observed.data))),
  ub.tvec = rep(Inf, times = length(as.numeric(observed.data))),
  std.error = FALSE,
  discrepancy = FALSE,
  user.ineq.constraint.function = NULL,
  opts.user = list(ftol_abs = 0, maxeval = 10000, xtol_rel = 1e-07, print_level = 0),
  zeroth.order = FALSE,
  method = c("constrained", "two_step"),
  control.nlminb = NULL
)
```

## Arguments

- observed.data:

  A numeric vector, matrix, or data frame containing the observed data.

- cgf:

  A CGF object corresponding to the distribution of the observed data.

- starting.theta:

  A numeric vector of starting values for model parameters.

- lb.theta:

  A numeric vector of lower bounds for `theta`. Defaults to `-Inf`.

- ub.theta:

  A numeric vector of upper bounds for `theta`. Defaults to `Inf`.

- starting.tvec:

  A numeric vector of starting values for the saddlepoint \\t\\.
  Defaults to `rep(0, length(observed.data))`.

- lb.tvec:

  Lower bounds for \\t\\. Defaults to `-Inf`.

- ub.tvec:

  Upper bounds for \\t\\. Defaults to `Inf`.

- std.error:

  Logical. If `TRUE`, compute standard errors of the MLEs. Defaults to
  `FALSE`.

- discrepancy:

  Logical. If `TRUE`, compute the discrepancy approximation. Defaults to
  `FALSE`.

- user.ineq.constraint.function:

  Optional user-defined inequality constraints on \\\theta\\. Must
  return a list with elements `constraints` and `jacobian`, with
  feasibility defined by `constraints <= 0`.

- opts.user:

  A named list of options for the nloptr optimizer. By default:
  `ftol_abs = 0`, `maxeval = 1e4`, `xtol_rel = 1e-7`, `print_level = 0`.
  The algorithm is fixed to `"NLOPT_LD_SLSQP"`.

- zeroth.order:

  Logical; if `TRUE` use the zeroth-order SPA objective. Defaults to
  `FALSE`.

- method:

  Optimisation strategy: `"constrained"` (default; backward compatible)
  or `"two_step"`.

- control.nlminb:

  Optional control list passed to `nlminb(..., control=...)` when
  `method="two_step"` and there are no user inequality constraints. If
  `NULL`, defaults are `eval.max=200`, `iter.max=150`, `rel.tol=1e-10`.

## Value

A list containing (at least):

- `MLEs.tvec, MLEs.theta`: the estimated saddlepoint \\\hat{t}\\ and
  parameter vector \\\hat{\theta}\\.

- `std.error`, `inverse.hessian`: if `std.error=TRUE`.

- `discrepancy`: if `discrepancy=TRUE`.

- `solution, status, message`: raw optimiser output (structure depends
  on optimiser).

## Details

- Observed Data: `observed.data` can be a numeric vector or a
  matrix/data frame. If it's a matrix/data frame, columns are treated as
  i.i.d. replicate blocks.

- Zeroth-order: If `zeroth.order = TRUE`, the objective omits the
  saddlepoint correction term and uses the zeroth-order approximation.

- Inequality constraints:

  - Constraints from `cgf$ineq_constraint(tvec,theta)` are enforced in
    `method="constrained"`.

  - In `method="two_step"`, `cgf$ineq_constraint` is enforced implicitly
    by the internal \\\hat t(\theta)\\ solver; candidates where no
    feasible \\\hat t(\theta)\\ can be found are penalized.

  - `user.ineq.constraint.function`, if supplied, is interpreted as a
    constraint on \\\theta\\ only. It must have signature
    `function(theta) -> list(constraints, jacobian)` with feasibility
    defined by `constraints <= 0` (NLOPT convention).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1); x <- rgamma(50, shape = 10, rate = 0.5)
res <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF,
                            starting.theta = c(1,1), std.error=TRUE)
res$MLEs.theta
res$std.error

## Two-step optimization:
res2 <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF,
                             starting.theta = c(1,1), method="two_step", std.error=TRUE)
res2$MLEs.theta
res2$std.error
} # }
```
