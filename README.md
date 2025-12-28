# saddlepoint

Composable cumulant generating function (CGF) objects and operators for saddlepoint approximations in R.

## What this package is for

`saddlepoint` provides:

- A CGF interface (K, K1, K2, …) for common distributions
- Operators to compose CGFs (i.i.d. replication, linear maps, sums of independent components, randomly stopped sums, …)
- Saddlepoint-based likelihood tools (e.g., MLE workflows) built on top of these CGFs

The pkgdown website contains the full reference and articles:
- https://godrick.github.io/saddlepoint/

## Installation

```r
# install.packages("pak")
pak::pak("godrick/saddlepoint")
```

### Quick start
```r
library(saddlepoint)

# Base scalar CGF example
lambda <- 3
PoissonCGF$K(tvec = 0.2, parameter_vector = lambda)
PoissonCGF$K1(tvec = 0.2, parameter_vector = lambda)
PoissonCGF$K2(tvec = 0.2, parameter_vector = lambda)

# CGF operators (examples)

set.seed(1)
cg_vec <- PoissonCGF
A <- matrix(c(1, 1), nrow = 1)
cg_sum <- linearlyMappedCGF(cg_vec, A)

lambda <- c(2, 3)
Y <- cg_sum$rsim(iidReps = 10000, parameter_vector = lambda)
mean(Y)  # ~ 5

# Randomly stopped sum (compound model)

set.seed(1)

theta <- c(lambdaN = 2, lambdaX = 3)

count_cgf   <- PoissonModelCGF(lambda = adaptor(indices = 1))
summand_cgf <- PoissonModelCGF(lambda = adaptor(indices = 2))

rss <- randomlyStoppedSumCGF(count_cgf, summand_cgf, block_size = 1)

Y <- rss$rsim(iidReps = 10000, parameter_vector = theta)
mean(Y)  # E[Y] = E[N]*E[X]
```

### CGF, derivatives and saddlepoint MLE
```R
# CGF at t = 0, shape = 10, rate = 0.5
GammaCGF$K(tvec = 0, parameter_vector = c(10, 0.5))

# First derivative (mean)
GammaCGF$K1(tvec = 0, parameter_vector = c(10, 0.5))

# Second derivative (variance)
GammaCGF$K2(tvec = 0, parameter_vector = c(10, 0.5))

# Sample data from Gamma distribution
set.seed(1); x = rgamma(50, shape = 10, rate = 0.5)

# MLE using saddlepoint likelihood
find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF, starting.theta = c(1,1))$MLEs.theta




