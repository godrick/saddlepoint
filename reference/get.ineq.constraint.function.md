# Create the inequality constraint function

This function wraps and returns an inequality constraint function of the
form `function(a) {...}`, where 'a' combines `tvec` and `theta`
arguments to a single vector, in that order.

## Usage

``` r
get.ineq.constraint.function(
  tvec,
  theta,
  cgf,
  user.ineq.constraint.function = NULL
)
```

## Arguments

- tvec:

  A numeric vector.

- theta:

  A numeric vector.

- cgf:

  An object of class 'CGF'.

- user.ineq.constraint.function:

  An optional additional inequality added by a user. Default is NULL.
  See TO DO list.

## Value

A function that takes a vector `a` as an argument. This function returns
either NULL or a list with 'constraints' and 'jacobian'.

## Details

This function constructs and integrates an inequality constraint based
on `cgf` with any user-defined inequality constraint function specified
in `user.ineq.constraint.function`. If `user.ineq.constraint.function`
is `NULL`, the function checks whether the `cgf` itself includes an
inherent inequality constraint.

## Examples

``` r
if (FALSE) { # \dontrun{
# TO DO: Add a working example
#  f <- get.ineq_constraint.function(tvec, theta, cgf)
#   f(c(tvec, theta)) # returns a list of the form list(constraints = , jacobian = ) or NULL
} # }
```
