---
title: saddlepoint
---

::: {.sp-hero}

# saddlepoint

`saddlepoint` provides a general framework for saddlepoint approximation (SPA) where distributions and model components are represented by cumulant generating functions (CGFs) and their derivatives. It supports saddlepoint likelihood evaluation and parameter estimation for composite and high-dimensional constructions.


Saddlepoint methods use the CGF
$
K(t) = \log \mathbb{E}\left[\exp(t X)\right]
$
to build accurate approximations to likelihoods, densities, and tail probabilities.

:::

## Key functionalities we provide

- CGF objects for common families (e.g., Poisson, Gamma, Gaussian, Binomial, Negative Binomial, Multinomial).
- Operators to build new CGFs by composition (e.g., i.i.d. replication, linear maps, sums, random sums, tilting operations, etc).
- Likelihood tools for saddlepoint-based inference, including workflows for maximum likelihood estimation.
- Error/discrepancy computations

## Quick start: saddlepoint MLE for a simple model

```r
library(saddlepoint)
set.seed(1)

y <- rnorm(60, mean = 2, sd = 1.5)

fit <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = NormalCGF,
  starting.theta = c(mu = 0, sigma = 1),
  lb.theta       = c(-Inf, 1e-8),
  ub.theta       = c( Inf, Inf),
  method         = "two_step",   
  std.error      = TRUE
)

fit$MLEs.theta
fit$std.error
```

For more detail (including constrained vs two-step fitting, a compositional aggregated-data example, and low-level likelihood evaluation), see: **[Saddlepoint likelihood inference](articles/mle.html)**.


## A simple example: MLE from aggregated data 

Consider this hypothetical scenario:

- Each day you have a random number of sessions $N \sim \mathrm{Poisson}(\lambda)$
- In each session there are $m$ opportunities (trials)
- Each opportunity succeeds with probability $p$
- You do not observe $N$. You only record the daily total successes
$$
  Y = \sum_{i=1}^N X_i, \quad X_i \sim \mathrm{Binomial}(m,p).
$$
For this simple model, the exact likelihood of $Y$ involves an infinite sum over the latent $N$:
$$
  \mathrm{Pr}(Y = y) = \sum_{n\ge 0} \mathrm{Pr}(N = n) \mathrm{Pr}(\mathrm{Binomial}(nm, p) = y)
$$

This is the kind of structure where SPA is convenient: we can build the CGF of $Y$ compositionally (using CGF building blocks) and fit $p$ without writing a custon pmf.


```r

library(saddlepoint)
set.seed(1)

B <- 60      # number of days
lambda <- 10  # mean sessions/day (fixed here for simplicity)
m <- 20      # opportunities per session
p_true <- 0.08

# Latent sessions per day 
N <- rpois(B, lambda)

# Observed: total successes per day
Y <- rbinom(B, size = N*m, prob = p_true)

# Count distribution for N (Poisson with fixed lambda)
count_cgf <- PoissonModelCGF(
  lambda  = adaptor(fixed_param = lambda)
)

# Summand distribution for X (Binomial with fixed m, unknown p)
# Here theta = p is 1-dimensional, so we map p via adaptor(indices = 1).
summand_cgf <- BinomialModelCGF(
  n = adaptor(fixed_param = m),
  p = adaptor(indices = 1)
)

# Randomly-stopped sum CGF for Y, replicated across B i.i.d. days
cgf <- randomlyStoppedSumCGF(
  count_cgf   = count_cgf,
  summand_cgf = summand_cgf,
  iidReps     = B
)

# Saddlepoint MLE for p
res_spa <- find.saddlepoint.MLE(
  observed.data = Y,
  cgf = cgf,
  starting.theta  = 0.01,
  lb.theta  = 1e-4,
  ub.theta  = 1 - 1e-4,
  std.error   = TRUE,
  discrepancy = TRUE
)
res_spa$MLEs.theta  
res_spa$std.error

# 'Exact' MLE via truncating the infinite sum (for comparison)
logSumExp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
log_pmf_Y <- function(y, p, lambda, m, n_max) {
  n_min <- if (y == 0) 0 else as.integer(ceiling(y / m))
  n <- n_min:n_max
  log_terms <- dpois(n, lambda = lambda, log = TRUE) +
    dbinom(y, size = n * m, prob = p, log = TRUE)
  logSumExp(log_terms)
}
n_max <- qpois(1 - 1e-12, lambda) + 50
nll_exact <- function(p) {
  if (!is.finite(p) || p <= 0 || p >= 1) return(Inf)
  -sum(vapply(
    Y, log_pmf_Y, numeric(1),
    p = p, lambda = lambda, m = m, n_max = n_max
  ))
}
opt_exact <- optimize(nll_exact, interval = c(1e-6, 1 - 1e-6))
p_hat_exact <- opt_exact$minimum
c(p_true = p_true,
  p_spa  = as.numeric(res_spa$MLEs.theta),
  p_exact_truncated = p_hat_exact)
  
res_spa$discrepancy # ~ p_hat_exact - res_spa$MLEs.theta
```
