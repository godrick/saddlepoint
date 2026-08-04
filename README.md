### `saddlepoint`

This package provides a saddlepoint approximation framework where distributions and model components are defined through cumulant generating functions (CGFs) and their derivatives. You can build new models by combining CGFs, and it supports likelihood-based inference via maximum likelihood estimation.

`saddlepoint` provides:

- A CGF interface (K, K1, K2, ...) for common distributions
- Operators to compose CGFs (i.i.d. replication, linear maps, sums of independent components, randomly stopped sums, ...)
- Saddlepoint-based likelihood tools (e.g., MLE workflows)

Examples and usage instructions are in the documentation and articles:
- https://godrick.github.io/saddlepoint/


### Installation

```r
# install.packages("devtools")
devtools::install_github("godrick/saddlepoint")
```

### Quick start
```r
library(saddlepoint)

# Base scalar CGF example
lambda <- 3
PoissonCGF$K(tvec = 0.2, parameter_vector = lambda)
PoissonCGF$K1(tvec = 0.2, parameter_vector = lambda)
PoissonCGF$K2(tvec = 0.2, parameter_vector = lambda)


# CGF objects
set.seed(1)
cg_vec <- PoissonCGF
A <- matrix(c(1, 1), nrow = 1)
cg_sum <- linearlyMappedCGF(cg_vec, A)

lambda <- c(10, 15)
Y <- cg_sum$rsim(n = 10000, vector_length = 1, parameter_vector = lambda)
mean(Y)  # ~ 25


# saddlepoint MLE
# sample data from Gamma distribution
set.seed(1); x = rgamma(50, shape = 10, rate = 0.5)
mle <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF, starting.theta = c(1,1))$MLEs.theta
mle 
```

### Citation

If you use `saddlepoint`, please cite the accompanying paper:

> Oketch, G., Fewster, R. M., and Goodman, J. (2026). *A general framework for computation and estimation using the saddlepoint approximation*. arXiv:2607.17464. <https://doi.org/10.48550/arXiv.2607.17464>

An automatically formatted citation and BibTeX entry are available in R with
`citation("saddlepoint")`.

