# Exponential Tilting of a CGF

Returns a new `CGF` corresponding to the exponential tilt of a base CGF.

If the base random vector has CGF \\K(t;\theta)\\, then the
exponentially tilted CGF is \$\$K\_{\text{tilt}}(t;\theta) = K(t +
h(\theta);\theta) - K(h(\theta);\theta).\$\$

## Usage

``` r
ExponentialTiltCGF(cgf, tilt, iidReps = "any", block_size = NULL, ...)
```

## Arguments

- cgf:

  A `CGF` object.

- tilt:

  A numeric vector, an `adaptor`, or a function mapping `theta -> h`.

  - If numeric, it is treated as a fixed tilt vector.

  - If an adaptor/function, it must return the tilt vector
    \\h(\theta)\\. The returned `h` must have length 1, or divide
    `length(tvec)` at evaluation time.

- iidReps:

  Either `"any"` or a positive integer.

- block_size:

  Optional block size for iid replication.

- ...:

  Passed to `createCGF` (advanced overrides).

## Value

A `CGF` object (tilted).

## Examples

``` r
## ------------------------------------------------------------
## Example 1: Poisson — tilt gives another Poisson (exact)
## ------------------------------------------------------------
lam <- 2
h   <- 0.4
cgT <- ExponentialTiltCGF(PoissonCGF, tilt = h, iidReps = 1)
tt <- c(-0.2, 0.1, 0.05)
stopifnot(all.equal(cgT$K(tt, lam), PoissonCGF$K(tt, lam * exp(h)), tol = 1e-12))

## ------------------------------------------------------------
## Example 2: Gamma — tilt shifts the rate (exact)
##   Gamma(a,b): K(tt) = -a log(1 - tt/b), domain tt < b
##   Tilt by h: Gamma(a, b-h), domain tt < (b-h)
## ------------------------------------------------------------
a <- 5; b <- 3; h <- 0.5
cgTg <- ExponentialTiltCGF(GammaCGF, tilt = h, iidReps = 1)
stopifnot(all.equal(cgTg$K(0.2, c(a,b)), GammaCGF$K(0.2, c(a, b - h)), tol = 1e-12))
# Check constraint: requires both tt+h < b and h < b
print(cgTg$ineq_constraint(0.2, c(a,b)))  # should be <= 0
#> [1] -2.3 -2.5

## ------------------------------------------------------------
## Example 3: No closed-form family - tilt a sum CGF
##   Z = Pois(lambda) + Binom(n,p)  (convolution, not a named family)
##   But we can still compute the exact likelihood numerically.
## ------------------------------------------------------------
set.seed(1)
n_fix <- 10
cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
cg_bin  <- BinomialModelCGF(n = adaptor(fixed_param = n_fix),
                            prob = adaptor(indices = 2), iidReps = 1)
cg_sum  <- sumOfIndependentCGF(list(cg_pois, cg_bin), iidReps = 1)

h <- 0.3
cg_sum_tilt <- ExponentialTiltCGF(cg_sum, tilt = h, iidReps = "any", block_size = 1)

# Simulate from the exact tilted model:
theta_base <- c(lambda = 2.5, p = 0.35)
lambda_t <- theta_base[1] * exp(h)
p_t      <- plogis(qlogis(theta_base[2]) + h)
B <- 50
y <- rpois(B, lambda_t) + rbinom(B, n_fix, p_t)

# Exact likelihood via convolution (finite sum):
dZ <- function(z, lambda, p) {
  k <- 0:n_fix
  sum(dbinom(k, n_fix, p) * dpois(z - k, lambda))
}
nll_exact <- function(theta) {
  lam <- theta[1]; p <- theta[2]
  lam_t <- lam * exp(h)
  p_t   <- plogis(qlogis(p) + h)
  -sum(log(vapply(y, dZ, numeric(1), lambda = lam_t, p = p_t)))
}

# SPA MLE using the tilted CGF
cgB <- iidReplicatesCGF(cg_sum_tilt, iidReps = B, block_size = 1)
fit_spa <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = cgB,
  starting.theta = c(1, 0.5),
  lb.theta       = c(1e-8, 1e-8),
  ub.theta       = c(Inf, 1-1e-8),
  method         = "two_step"
)

# Exact MLE (numeric) for comparison:
fit_exact <- nlminb(
  start = c(1, 0.5),
  objective = nll_exact,
  lower = c(1e-8, 1e-8),
  upper = c(Inf, 1 - 1e-8)
)

cat("True theta_base:", theta_base, "\n")
#> True theta_base: 2.5 0.35 
cat("SPA  theta_hat :", fit_spa$MLEs.theta, "\n")
#> SPA  theta_hat : 1.440063 0.5073478 
cat("Exact theta_hat:", fit_exact$par, "\n")
#> Exact theta_hat: 1.342831 0.5208833 

## ------------------------------------------------------------
## Example 4: theta-dependent tilt + SPA MLE vs exact MLE
##   Z = Pois(lambda) + Binom(n,p), then exponential-tilt by h(theta)
##   where h(theta) = log(lambda).
##
##   Under tilt by h, the components tilt as:
##     Pois(lambda)    -> Pois(lambda * exp(h))
##     Binom(n,p)      -> Binom(n, plogis(qlogis(p) + h))
##
##   Here h depends on theta, so the tilted likelihood is still a
##   well-defined model in theta, but it is NOT “just a fixed tilted
##   family” anymore.
## ------------------------------------------------------------
# \donttest{
set.seed(2)
B     <- 60
n_fix <- 10

## Base CGF for one observation Z = X + Y
cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
cg_bin  <- BinomialModelCGF(
  n    = adaptor(fixed_param = n_fix),
  prob = adaptor(indices = 2),
  iidReps = 1
)
cg_sum <- sumOfIndependentCGF(list(cg_pois, cg_bin), iidReps = 1)

## Theta-dependent tilt:
## h(theta) = log(lambda).  (Requires lambda>0; we enforce via lb.theta.)
tilt_theta <- function(theta) log(theta[1])

## tilted CGF for B i.i.d. observations (scalar blocks)
cg_tilt_B <- ExponentialTiltCGF(
  cgf        = cg_sum,
  tilt       = tilt_theta,
  iidReps    = B,
  block_size = 1
)

## Simulate data from the exact tilted model at theta_true:
theta_true <- c(lambda = 2.5, p = 0.35)
h_true     <- tilt_theta(theta_true)

lambda_t <- theta_true[1] * exp(h_true)                 # = lambda^2
p_t      <- plogis(qlogis(theta_true[2]) + h_true)

y <- rpois(B, lambda_t) + rbinom(B, n_fix, p_t)

## SPA MLE under the tilted CGF
fit_spa <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = cg_tilt_B,
  starting.theta = c(1.5, 0.5),
  lb.theta       = c(1e-8, 1e-8),
  ub.theta       = c(Inf, 1 - 1e-8),
  method         = "two_step"
)

## finite convolution for comparison
dZ <- function(z, lambda, p) {
  k <- 0:n_fix
  sum(dbinom(k, n_fix, p) * dpois(z - k, lambda))
}

nll_exact <- function(theta) {
  h       <- tilt_theta(theta)
  lam_t   <- theta[1] * exp(h)
  p_tilt  <- plogis(qlogis(theta[2]) + h)
  likvec  <- vapply(y, dZ, numeric(1), lambda = lam_t, p = p_tilt)
  -sum(log(pmax(likvec, .Machine$double.xmin)))
}

fit_exact <- nlminb(
  start     = c(1.5, 0.5),
  objective = nll_exact,
  lower     = c(1e-8, 1e-8),
  upper     = c(Inf, 1 - 1e-8)
)

cat("theta_true        =", paste(round(theta_true, 4), collapse=" "), "\n")
#> theta_true        = 2.5 0.35 
cat("theta_hat_spa     =", paste(round(fit_spa$MLEs.theta, 4), collapse=" "), "\n")
#> theta_hat_spa     = 2.9718 0.1516 
cat("theta_hat_exact   =", paste(round(fit_exact$par, 4), collapse=" "), "\n")
#> theta_hat_exact   = 2.9071 0.1771 
# }
```
