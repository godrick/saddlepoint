# Create an adapted CGF object using an adaptor function

Constructs a new `CGF` object by adapting the parameter vector using a
user-supplied `adaptor` before invoking the original `CGF`'s methods.

## Usage

``` r
adaptCGF(cgf, adaptor, ...)
```

## Arguments

- cgf:

  A `CGF` object to be adapted.

- adaptor:

  An `adaptor` or a function with signature
  `function(theta) -> adapted_param`.

- ...:

  Additional named arguments passed to
  [`createCGF`](https://godrick.github.io/saddlepoint/reference/createCGF.md),
  the `CGF` object creation function (rarely needed).

## Value

A new `CGF` object.

## Details

This function is useful when you have a `CGF` that expects its parameter
vector to be in a certain format, but your high-level model provides a
different parameter structure. By specifying a `adaptor`, you can
dynamically translate your model parameters to those that `cgf`
requires.

## Examples

``` r
if (FALSE) { # \dontrun{
## Example: Suppose you have a sum of two Poisson r.v.s, each with
##   a different lambda, but your model uses a single parameter 'theta'
##   from which the two lambdas are derived (lambda1 = theta, lambda2 = 2*theta).

## Scenario:
## Y1 ~ Poisson(theta)
## Y2 ~ Poisson(2*theta)
## Y = Y1 + Y2 ~ Poisson(3*theta)
## Model Parameters: theta
## CGF for Y expects: distribution_params = c(lambda1, lambda2) = c(theta, 2*theta)


# First, the individual Poisson CGFs with separate lambda parameters
# These will expect a vector of length 2 for the two lambdas
K_Y1 <- PoissonModelCGF(lambda = adaptor(indices = 1))   # For Y1: lambda1 => the first parameter
K_Y2 <- PoissonModelCGF(lambda = adaptor(indices = 2))   # For Y2: lambda2 => the second parameter

# sum_cgf is a CGF that expects 2 parameters for the two Poisson variables
sum_cgf <- sumOfIndependentCGF(cgf_list = list(K_Y1, K_Y2))

# adaptor that converts 'theta' to c(lambda1, lambda2)
mapThetaToDistParams <- function(theta) c(theta, 2*theta)

# Now adapt sum_cgf so it only needs a single 'theta':
adapted_cgf <- adaptCGF(cgf = sum_cgf, adaptor = mapThetaToDistParams)
theta0 <- 5
adapted_cgf$K1(0, theta0)
# OR sum_cgf$K1(0, c(theta0, 2*theta0))
# OR PoissonCGF$K1(0, 3*theta0)
} # }
```
