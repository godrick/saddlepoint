---
title: saddlepoint
---

<div class="sp-hero">
<div class="container">

# saddlepoint

Composable cumulant generating function (CGF) objects for saddlepoint approximations in R.

Saddlepoint methods use the CGF
\[
K(t) = \log \mathbb{E}\left[\exp(t^\top X)\right]
\]
to build accurate approximations to likelihoods, densities, and tail probabilities.

</div>
</div>

## What this package provides

- **CGF objects** for common families (e.g., Poisson, Gamma, Gaussian, Binomial, Negative Binomial, Multinomial).
- **Operators** to build new CGFs by composition (e.g., i.i.d. replication, linear maps, sums, random sums, Esscher tilts).
- **Likelihood tools** for saddlepoint-based inference, including workflows for maximum likelihood estimation.

## A small example

```r
library(saddlepoint)
set.seed(1)

# Sum of two independent Poisson components:
# X = (X1, X2) with lambda = (2, 3), and Y = X1 + X2.
A <- matrix(c(1, 1), nrow = 1)
cg_sum <- linearlyMappedCGF(PoissonCGF, A)

lambda <- c(2, 3)
y <- cg_sum$rsim(iidReps = 10000, parameter_vector = lambda)

c(mean = mean(y), var = var(y))
