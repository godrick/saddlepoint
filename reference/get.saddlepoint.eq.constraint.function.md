# Create the saddlepoint equality constraint function

This function creates and returns a function of the form
`function(a) {...}`, where 'a' combines `tvec` and `theta` arguments to
a single vector, in that order.

## Usage

``` r
get.saddlepoint.eq.constraint.function(tvec, theta, observed.data, cgf)
```

## Arguments

- tvec:

  A numeric vector.

- theta:

  A numeric vector.

- observed.data:

  A numeric vector of observations with the same length as `tvec`.

- cgf:

  An object of class 'CGF'.

## Value

A function that accepts a vector 'a' as an argument. When
`a = c(tvec, theta)` is passed to this function, it generates a list
containing 'constraints' and 'jacobian'. 'constraints' are computed as
\\K'(t;\theta) - y\\, and 'jacobian' represents the gradient of these
constraints with respect to both `tvec` and `theta`.

## Details

In this set up, the expression of the saddlepoint equation is defined by
\\K'(t;\theta) - y = 0\\ for `observed.data = y`. The returned function
is designed to compute both \\K'(t;\theta) - y\\ and its gradient with
respect to both `tvec` and `theta`.
