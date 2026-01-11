# Shifted CGF (deterministic translation)

Returns the CGF of \\Y = X + b(\theta)\\, where \\X\\ has CGF `cgf` and
\\b(\theta)\\ is a deterministic shift (constant or theta-dependent).

This operation does not change the dimension of the random vector and
therefore does not impose any new i.i.d. replication semantics. If you
want a replicated interpretation of `tvec`, wrap the result with
[`iidReplicatesCGF`](https://godrick.github.io/saddlepoint/reference/iidReplicatesCGF.md).

## Usage

``` r
shiftedCGF(cgf, shift, ...)
```

## Arguments

- cgf:

  A `CGF` object.

- shift:

  A numeric scalar/vector, an `adaptor`, or a function `shift(theta)`
  returning the shift vector. If its length is 1, it is recycled; if its
  length divides `length(tvec)`, it is recycled blockwise.

- ...:

  Passed to
  [`createCGF()`](https://godrick.github.io/saddlepoint/reference/createCGF.md)
  (advanced).

## Value

A `CGF` object.

## Examples

``` r
## Example 1: Poisson shifted by a constant (exact relationship)
lam <- 2
b   <- 3
cgS <- shiftedCGF(PoissonCGF, shift = b)
tt   <- 0.2
stopifnot(all.equal(cgS$K(tt, lam), PoissonCGF$K(tt, lam) + b*tt, tol = 1e-12))
stopifnot(all.equal(cgS$K1(tt, lam), PoissonCGF$K1(tt, lam) + b, tol = 1e-12))

## Example 2: Multivariate Normal with theta-dependent shift + MLE
if (FALSE) { # \dontrun{
set.seed(1)
d <- 2
B <- 200

Sigma <- matrix(c(1.0, 0.3,
                  0.3, 2.0), nrow = d, byrow = TRUE)
mu0 <- c(-1, 0.2)

# Base MVN with fixed mean/covariance (use fixed_param to keep AD context)
cg_base <- MultivariateNormalModelCGF(
  mu    = adaptor(fixed_param = mu0),
  sigma = adaptor(fixed_param = Sigma),
  iidReps = "any"
)

# Shift b(theta) in R^d
b_fun <- function(theta) c(theta[1], exp(theta[2]))
cg <- shiftedCGF(cg_base, shift = b_fun)

# Simulate Y = X + b(theta_true)
theta_true <- c(0.7, log(2))     # b(theta_true) = (0.7, 2)
mu_true <- mu0 + b_fun(theta_true)
Y <- cg$rsim(n = B, vector_length = d, parameter_vector = theta_true)

fit <- find.saddlepoint.MLE(
  observed.data  = Y,            # columns treated as i.i.d. blocks
  cgf            = cg,
  starting.theta = c(0, 0),
  method         = "two_step"
)
fit$MLEs.theta

# Closed-form MLE (known Sigma): match the sample mean
m <- rowMeans(Y)
c(m[1] - mu0[1], log(m[2] - mu0[2]))
} # }

```
