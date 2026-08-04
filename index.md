# `saddlepoint`

`saddlepoint` provides a general framework for saddlepoint approximation
(SPA) where distributions and model components are represented by
cumulant generating functions (CGFs) and their derivatives. It supports
saddlepoint likelihood evaluation and parameter estimation for composite
and high-dimensional constructions.

Saddlepoint methods use the CGF
$`K(t) = \log \mathbb{E}\left[\exp(t X)\right]`$ to build accurate
approximations to likelihoods, densities, and tail probabilities.

## Key functionalities we provide

- CGF objects for common families (e.g., Poisson, Gamma, Gaussian,
  Binomial, Negative Binomial, Multinomial).
- Operators to build new CGFs by composition (e.g., i.i.d. replication,
  linear maps, sums, random sums, tilting operations, etc).
- Likelihood tools for saddlepoint-based inference, including workflows
  for maximum likelihood estimation.
- Error/discrepancy computations

## Quick start: saddlepoint MLE for a simple model

``` r


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

## Citation

If you use `saddlepoint`, please cite the accompanying paper:

> Oketch, G., Fewster, R. M., and Goodman, J. (2026). *A general
> framework for computation and estimation using the saddlepoint
> approximation*. arXiv:2607.17464.
> <https://doi.org/10.48550/arXiv.2607.17464>

Use `citation("saddlepoint")` for an automatically formatted citation
and BibTeX entry.
