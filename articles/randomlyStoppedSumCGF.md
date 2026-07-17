# Randomly-stopped sums

This vignette explains how to build and use a CGF for a randomly-stopped
sum with
[`randomlyStoppedSumCGF()`](https://godrick.github.io/saddlepoint/reference/randomlyStoppedSumCGF.md),
and (most importantly) how to avoid the most common replication/blocking
mistakes.

## Model

A **randomly-stopped sum** has the form

``` math
Y = \sum_{i=1}^{N} X_i,
```

where

- $`N`$ is a scalar non-negative integer count random variable,
- $`X_1,X_2,\dots`$ are i.i.d. copies of a (possibly vector-valued)
  random vector $`X`$,
- $`N`$ is independent of all $`X_i`$.

In CGF terms, the key identity is the composition

``` math
K_Y(t;\theta) = K_N\bigl( K_X(t;\theta);\theta \bigr),
```

where $`K_X(t;\theta)`$ is a scalar even when $`t`$ is a vector.

You provide:

- `count_cgf`: a CGF for $`N`$ (scalar),
- `summand_cgf`: a CGF for $`X`$ (scalar or vector),
- a shared parameter vector `theta` (both CGFs read from the same
  `theta`).

``` r


library(saddlepoint)
```

## The one thing to get right: `block_size` and `iidReps`

Many CGFs in this package are vectorized and can interpret a long `tvec`
in more than one way. For RSS models that ambiguity can silently change
the statistical model (scalar vs vector, “one observation” vs “many
i.i.d. observations”).

[`randomlyStoppedSumCGF()`](https://godrick.github.io/saddlepoint/reference/randomlyStoppedSumCGF.md)
forces you to disambiguate:

- `block_size` = dimension of one observation $`Y_b`$ (equals dim of the
  summand $`X`$).
- `iidReps` = number of i.i.d. observations $`B`$ in your dataset.

### Recommended patterns

Scalar RSS data (a numeric vector `y` of length `B`):

- safest: specify `block_size = 1` and let `iidReps = "any"` (defaulted
  internally).

``` r


cg_rss <- randomlyStoppedSumCGF(count_cgf   = PoissonCGF,
                                summand_cgf = BinomialModelCGF(n = adaptor(fixed_param = 1),
                                                               prob = adaptor(fixed_param = 0.3)), 
                                block_size  = 1)
```

Vector RSS data (each observation is length `d`, stored as a `d x B`
matrix):

- set `block_size = d`,
- either set `iidReps = "any"` (and let it infer `B` from
  `length(tvec) / d`) or set `iidReps = B`.

### What happens if you get it wrong?

If you forget `block_size` / `iidReps`, a “long” `tvec` might be treated
as:

- a higher-dimensional summand $`X`$, which changes the model for $`Y`$,
  or
- multiple i.i.d. blocks in some child CGF (again changing $`K_X(t)`$,
  hence changing $`K_Y(t)`$).

That is why RSS requires at least one of `block_size` or `iidReps`.

## Example 1: Poisson count + Bernoulli summand = Poisson thinning (closed form)

Let

- $`N \sim \text{Poisson}(\lambda)`$,
- $`X \sim \text{Bernoulli}(p)`$, i.i.d.,
- $`Y = \sum_{i=1}^N X_i`$.

Then $`Y \sim \text{Poisson}(\lambda p)`$ (Poisson thinning), so the RSS
CGF should match a Poisson CGF with mean $`\lambda p`$.

``` r


lambda <- 4
p      <- 0.3
theta  <- c(lambda = lambda, p = p)

cg_N <- PoissonModelCGF(lambda = adaptor(indices = 1))
cg_X <- BinomialModelCGF(n = adaptor(fixed_param = 1), prob = adaptor(indices = 2))

cg_rss <- randomlyStoppedSumCGF(count_cgf   = cg_N,
                               summand_cgf = cg_X,
                               block_size  = 1) # scalar Y

tt <- c(-0.2, 0.1, 0.05)

# RSS CGF:
K_rss <- cg_rss$K(tt, theta)

# Closed-form equivalent:
K_closed <- PoissonCGF$K(tt, lambda * p)

cbind(K_rss, K_closed)
#>            K_rss    K_closed
#> [1,] -0.02979268 -0.02979268
stopifnot(all.equal(K_rss, K_closed, tol = 1e-12))
```

## Example 2: Poisson count + categorical summand = Poisson splitting (vector, closed form)

Let

- $`N \sim \text{Poisson}(\lambda)`$,
- $`X`$ is a ‘one-hot’ vector in $`\{e_1,\dots,e_d\}`$ with
  probabilities $`\pi`$,
- $`Y = \sum_{i=1}^N X_i`$ is a vector of category counts.

Then the components of $`Y`$ are independent Poisson with means
$`\lambda \pi_j`$.

``` r


set.seed(1)
d <- 3
lambda <- 10
pi_ <- c(0.2, 0.5, 0.3)
theta <- c(lambda = lambda, pi_)

cg_N <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)

# summand: Multinomial with n=1 gives a categorical one-hot vector
cg_X <- MultinomialModelCGF(
  n = adaptor(fixed_param = 1),
  prob_vec = adaptor(indices = 2:(d+1)),
  iidReps = 1
)

cg_rss_vec <- randomlyStoppedSumCGF(cg_N, cg_X, block_size = d)

# Check K1 at t=0 equals E[Y] = lambda * pi_
cg_rss_vec$K1(rep(0, d), theta)
#> [1] 2 5 3
lambda * pi_
#> [1] 2 5 3
```

## Example 3: Compound Poisson–Gamma

A classic ‘not-nicely-closed-form’ RSS is the compound Poisson–Gamma:

- $`N \sim \text{Poisson}(\lambda)`$,
- $`X \sim \text{Gamma}(\alpha,\beta)`$ i.i.d. severities (shape
  $`\alpha`$, rate $`\beta`$),
- $`Y = \sum_{i=1}^N X_i`$.

The distribution of $`Y`$ has a point mass at 0 and a continuous density
for $`y>0`$. The exact density is an infinite series
``` math
f_Y(y) = e^{-\lambda}\sum_{k=1}^{\infty}\frac{\lambda^k}{k!} f_{\Gamma}(y;\,k\alpha,\beta),
```
so in practice “exact” likelihood evaluation typically needs truncation
or special-purpose methods.

This is a good use case for saddlepoint likelihoods.

### Build the RSS CGF

``` r


set.seed(123)

alpha_fix <- 2
theta_true <- c(lambda = 10, rate = 3)

cg_N <- PoissonModelCGF(lambda = adaptor(indices = 1) ) # optional # iidReps = 1
cg_X <- GammaModelCGF(shape = adaptor(fixed_param = alpha_fix),
                      rate  = adaptor(indices = 2)) # optional # iidReps = 1

cg_rss <- randomlyStoppedSumCGF(cg_N, cg_X, block_size = 1)
```

### Simulate data

``` r


B <- 100
Y <- as.numeric(cg_rss$rsim(n = B, vector_length = 1, parameter_vector = theta_true))
summary(Y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.099   4.856   5.988   6.183   7.644  11.917
```

### Fit via saddlepoint likelihood

``` r


fit_spa <- find.saddlepoint.MLE(
  observed.data  = Y,
  cgf            = cg_rss,  # iidReps="any" because we supplied block_size=1
  starting.theta = c(2, 1),
  lb.theta       = c(1e-3, 1e-6),
  method         = "constrained",  # usually safest when there are domain constraints
  std.error      = TRUE,
  discrepancy = TRUE
)

fit_spa$MLEs.theta
#> [1] 10.208868  3.302203
fit_spa$std.error
#> [1] 1.4671329 0.4689042
```

### Compare against an exact likelihood

The following uses the infinite-series density truncated at `kmax`.

``` r


d_comp_pois_gamma <- function(y, lambda, alpha, rate, kmax = 200) {
  if (y == 0) return(exp(-lambda))
  ks <- 1:kmax
  weights <- exp(-lambda) * lambda^ks / factorial(ks)
  dens <- dgamma(y, shape = ks * alpha, rate = rate)
  sum(weights * dens)
}

nll_trunc <- function(theta) {
  lam <- theta[1]; rate <- theta[2]
  -sum(log(vapply(Y, d_comp_pois_gamma, numeric(1),
                  lambda = lam, alpha = alpha_fix, rate = rate, kmax = 200)))
}

fit_trunc <- nlminb(
  start = c(1, 1),
  objective = nll_trunc,
  lower = c(1e-6, 1e-6)
)

rbind(true = theta_true, spa = fit_spa$MLEs.theta, true_mle = fit_trunc$par)
fit_spa$discrepancy                # approximated difference
fit_trunc$par - fit_spa$MLEs.theta # exact difference
```

## Practical checklist

When using
[`randomlyStoppedSumCGF()`](https://godrick.github.io/saddlepoint/reference/randomlyStoppedSumCGF.md):

1.  Always decide what one observation is.
    - If one observation is scalar, set `block_size = 1`.  
    - If one observation is a vector of length `d`, set
      `block_size = d`.
2.  Decide how you store the data.
    - scalar data: `observed.data` as a vector of length `B`  
    - vector data: `observed.data` as a `d x B` matrix (columns = i.i.d.
      blocks)
3.  If constraints exist, prefer `method="constrained"` in
    [`find.saddlepoint.MLE()`](https://godrick.github.io/saddlepoint/reference/find.saddlepoint.MLE.md).
