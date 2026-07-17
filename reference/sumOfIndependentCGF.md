# CGF Object for the sum of independent random variables

Constructs a new CGF object representing the sum of independent random
vectors: \\Y = Y^{(1)} + \cdots + Y^{(L)}\\ where each summand has its
own CGF and the summands are independent.

## Usage

``` r
sumOfIndependentCGF(cgf_list, iidReps = NULL, block_size = NULL, ...)
```

## Arguments

- cgf_list:

  A non-empty list of CGF objects (each inherits from class `"CGF"`).

- iidReps:

  Optional. Either `NULL`, `"any"`, or a positive integer.

- block_size:

  Either `NULL` or a positive integer describing the block size for iid
  replication.

- ...:

  Additional named arguments passed to
  [`createCGF`](https://godrick.github.io/saddlepoint/reference/createCGF.md)
  (rare).

## Value

A `CGF` object.

## Details

**Replication (iidReps / block_size):**

- If both `iidReps` and `block_size` are `NULL`, no replication is
  applied.

- If `block_size` is provided but `iidReps` is `NULL`, we set
  `iidReps = "any"` and infer the number of blocks from
  `length(tvec) / block_size` at evaluation time.

- If `iidReps = "any"`, then `block_size` must be provided.

- If `iidReps` is a positive integer, `block_size` may be `NULL`, though
  providing it is encouraged.

Note: `iidReps`/`block_size` describe i.i.d. replication of the sum
\\Y\\, not the length of `cgf_list`.

## Examples

``` r
## -----------------------------
## Sum of independent Poisson variables
## -----------------------------
## Here we build Y = Y1 + Y2 with Y1 ~ Pois(lambda1), Y2 ~ Pois(lambda2),
## independent. The mean is lambda1 + lambda2, i.e. K1(0) should equal that.
lambda1 <- 2
lambda2 <- 5

cg1 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda1), iidReps = 1)
cg2 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda2), iidReps = 1)

## one observation (block_size = 1)
cg_sum <- sumOfIndependentCGF(list(cg1, cg2), iidReps = 1, block_size = 1)

theta_dummy <- 0
cg_sum$K1(tvec = 0, parameter_vector = theta_dummy)   # should be 7
#> [1] 7


## ------------------------------------------------------------
## Sum of independent Gammas with different rates
## (demonstrates concatenated inequality constraints)
## ------------------------------------------------------------
## If X ~ Gamma(shape=a, rate=b), its CGF is only valid for tvec < b.
## For Y = X1 + X2 with rates b1 and b2, the domain is t < min(b1, b2),
## and the implementation is such that the constraint vector is concatenated:
##   g(t,theta) = c(t - b1, t - b2)  <= 0
##
# \donttest{
set.seed(123)
B <- 50

## True shapes (fixed) and rates (unknown; theta = (b1, b2))
a1_true <- 10
a2_true <- 3
b1_true <- 7
b2_true <- 1
theta_true <- c(b1_true, b2_true)

## Two GammaModelCGFs, each uses one component of theta as its rate
cg_g1 <- GammaModelCGF(
  shape  = adaptor(fixed_param = a1_true),
  rate   = function(th) th[1],
  iidReps = 1
)
cg_g2 <- GammaModelCGF(
  shape  = adaptor(fixed_param = a2_true),
  rate   = function(th) th[2],
  iidReps = 1
)

## Sum for ONE observation
cg_sum_one <- sumOfIndependentCGF(list(cg_g1, cg_g2), iidReps = 1, block_size = 1)

## Inequality constraints are concatenated across summands:
## For scalar t: g(t,theta) = c(t-b1, t-b2) must be <= 0.
cg_sum_one$ineq_constraint(tvec = 0.9, param = theta_true)  # ~ c(-6.1, -0.1) (feasible)
#> [1] -6.1 -0.1
cg_sum_one$ineq_constraint(tvec = 1.2, param = theta_true)  # ~ c(-5.8, +0.2) (violates 2nd)
#> [1] -5.8  0.2

## B i.i.d. replicates of the sum
cg_sum_B <- sumOfIndependentCGF(list(cg_g1, cg_g2), iidReps = B, block_size = 1)

## Simulate data y_i = x1_i + x2_i
y <- rgamma(B, shape = a1_true, rate = b1_true) +
     rgamma(B, shape = a2_true, rate = b2_true)

## NOTE / recommendation:
## When cgf$ineq_constraint is non-empty (domain-constrained CGFs),
## method="constrained" is currently the most robust choice because it enforces
## the CGF domain constraints directly during optimisation. The "two_step" method
## WILL be slower for constrained CGFs in the current implementation.

fit_const <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = cg_sum_B,
  starting.theta = c(1.2, 0.5),
  lb.theta       = c(1e-4, 2e-5),
  method         = "constrained"
)

fit_two <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = cg_sum_B,
  starting.theta = c(1.2, 0.5),
  lb.theta       = c(1e-4, 2e-5),
  method         = "two_step"
)
#> Warning: method='two_step' was requested, but this CGF has inequality constraints on tvec. In the current implementation, the two-step approach will be extremely slow. For constrained CGFs, consider method='constrained' for speed/robustness. Performance for this case may improve in future versions.
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced

## Optional quick checks: saddlepoint residual and feasibility
max(abs(cg_sum_B$K1(fit_const$MLEs.tvec, fit_const$MLEs.theta) - y))
#> [1] 3.552714e-15
max(cg_sum_B$ineq_constraint(fit_const$MLEs.tvec, fit_const$MLEs.theta))  # should be <= 0
#> [1] -0.4579342

cat("true theta:", theta_true, "\n")
#> true theta: 7 1 
cat("constrained:", round(fit_const$MLEs.theta, 4), "\n")
#> constrained: 4.9247 1.3149 
cat("two_step   :", round(fit_two$MLEs.theta, 4), "\n")
#> two_step   : 4.9247 1.3149 
# }

```
