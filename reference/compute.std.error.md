# Compute standard error and inverse Hessian

This function calculates the standard error and the inverse of the
Hessian matrix for a model incorporating both the saddlepoint likelihood
and an optional non-saddlepoint negative log-likelihood component.

## Usage

``` r
compute.std.error(
  observed.data,
  estimated.tvec,
  estimated.theta,
  cgf,
  zeroth.order = FALSE,
  non.saddlepoint.negll.function = NULL
)
```

## Arguments

- observed.data:

  A numeric vector of observations used in the saddlepoint likelihood
  calculations.

- estimated.tvec:

  A numeric vector of MLE values for the saddlepoint values `tvec`.

- estimated.theta:

  A numeric vector of MLE values for the model parameters.

- cgf:

  A CGF object used in the saddlepoint likelihood computation.

- zeroth.order:

  A logical value indicating whether the zeroth-order saddlepoint
  likelihood is used. Default is FALSE.

- non.saddlepoint.negll.function:

  An optional function specifying a non-saddlepoint negative
  log-likelihood. If your model used for estimation incorporates a
  likelihood component that is not based on the saddlepoint
  approximation, this function must be provided. See details for more
  information.

## Value

A list containing:

- std.error: Standard error computed from the diagonal of the inverse
  Hessian.

- inverse.hessian: Inverse of the Hessian matrix.

## Details

The `non.saddlepoint.negll.function`, if provided, should have the form:
\$\$function(theta) \\ return(list(objective = ..., hessian = ...))
\\\$\$ It adds flexibility by allowing incorporation of likelihood
components not based on the saddlepoint likelihood. It is important when
the analysis combines the saddlepoint likelihood with an external
likelihood. If not provided, the function defaults to using solely the
saddlepoint likelihood.
