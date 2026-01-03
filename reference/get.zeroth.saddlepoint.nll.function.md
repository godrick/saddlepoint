# Create the zeroth-order saddlepoint negative log-likelihood function

This function creates and returns a function of the form
`function(a) {...}`, where 'a' combines `tvec` and `theta` arguments to
a single vector, in that order.

## Usage

``` r
get.zeroth.saddlepoint.nll.function(tvec, theta, cgf)
```

## Arguments

- tvec:

  A numeric vector.

- theta:

  A numeric vector.

- cgf:

  An object of class 'CGF'.

## Value

A function that takes a vector 'a' as an argument. When
`a = c(tvec, theta)` is passed to the returned function, it yields a
list in the format: `list(objective = , gradient = )`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # TO DO: write a working example
} # }
```
