# Create a CGF object from user-defined functions

This creates an object of type CGF using user-supplied functions. You
supply the five essential methods (`K`, `K1`, `K2`, `K3operator`,
`K4operator`) plus any optional overrides (e.g., `tilting_exponent` or
`neg_ll`), and it returns a `CGF` instance.

## Usage

``` r
createCGF(
  K,
  K1,
  K2,
  K3operator,
  K4operator,
  ineq_constraint = NULL,
  analytic_tvec_hat_func = NULL,
  op_name = "UnnamedOperation",
  tilting_exponent = NULL,
  neg_ll = NULL,
  func_T = NULL,
  K2_solve = NULL,
  logdetK2 = NULL,
  rsim = NULL,
  K4operatorAABB = NULL,
  K3K3operatorAABBCC = NULL,
  K3K3operatorABCABC = NULL,
  K4operatorAABB_factored = NULL,
  K3K3operatorAABBCC_factored = NULL,
  K3K3operatorABCABC_factored = NULL,
  K2operator = NULL,
  K2operatorAK2AT = NULL,
  ...
)
```

## Arguments

- K:

  A function `K(tvec, parameter_vector) -> numeric scalar`.

- K1:

  A function `K1(tvec, parameter_vector) -> numeric vector`.

- K2:

  A function `K2(tvec, parameter_vector) -> numeric matrix`.

- K3operator:

  A function implementing the third-order operator.

- K4operator:

  A function implementing the fourth-order operator.

- ineq_constraint:

  Optional function for inequality constraints.

- analytic_tvec_hat_func:

  Optional function for an analytic solution of the saddlepoint
  equation. If provided, call it via `cgf$analytic_tvec_hat(x, param)`.

- op_name:

  A descriptive label for the CGF object/operation. Default is
  "UnnamedOperation".

- tilting_exponent:

  (optional) Overriding function for the tilting exponent.

- neg_ll:

  (optional) Overriding function for the negative log-likelihood.

- func_T:

  (optional) Overriding function for the first-order correction term.

- K2_solve:

  Optional computational helper. A function of the form
  `function(tvec, parameter_vector, rhs)` returning the solution to
  `K2(tvec, parameter_vector) %*% x = rhs`. `rhs` may be a vector or a
  matrix (solve column-wise).

- logdetK2:

  Optional computational helper. A function of the form
  `function(tvec, parameter_vector)` returning
  `log(det(K2(tvec, parameter_vector)))` Useful for wrapper CGFs that
  can compute this without materializing the full `K2`.

- rsim:

  Optional simulation method. A function of the form
  `function(iidReps, parameter_vector)` returning a numeric vector or
  matrix of simulated observations. If supplied, the resulting CGF
  exposes `$rsim()` and `$has_simulate()`.

- K4operatorAABB_factored, K3K3operatorAABBCC_factored,
  K3K3operatorABCABC_factored:

  (optional) Overriding factored-operator methods.

- K2operator, K2operatorAK2AT, K4operatorAABB, K3K3operatorAABBCC,
  K3K3operatorABCABC:

  (optional) Overriding operator methods.

- ...:

  Additional named methods or overrides.

## Value

An object of class `CGF`.
