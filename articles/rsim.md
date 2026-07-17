# Random simulation with CGF objects: rsim

``` r


library(saddlepoint)
```

This vignette is about `CGF$rsim()`, which simulates draws from the
distribution described by a `CGF` object.

## `rsim()` arguments

All `CGF$rsim()` methods follow the same return convention:

- `n` is the number of independent draws to generate.
- `vector_length` is the dimension of one draw (each draw is a column).
- The return value is a matrix of shape `vector_length x n`.
- With `flatten = TRUE`, the result is flattened to a numeric vector
  (column-major order).

``` r


set.seed(1)

lambda <- 2
X <- PoissonCGF$rsim(n = 5, vector_length = 1, parameter_vector = lambda)
X
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    1    2    4    1
dim(X)
#> [1] 1 5
```

If you prefer a plain vector:

``` r


PoissonCGF$rsim(n = 5, vector_length = 1, parameter_vector = lambda, flatten = TRUE)
#> [1] 4 4 2 2 0
```

## Vector-valued simulation

CGF objects support vector-valued random variables. For example, a
vector of independent Poisson variables can be represented by passing a
vector `parameter_vector` (one `lambda` per component) and setting
`vector_length` to match.

``` r


set.seed(1)

lambda_vec <- c(2, 5, 10)
Xv <- PoissonCGF$rsim(n = 4, vector_length = length(lambda_vec), parameter_vector = lambda_vec)
Xv
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    2    1
#> [2,]    4    8    2    7
#> [3,]   10   15    7    9
dim(Xv)
#> [1] 3 4
```

### Linear maps + replicated blocks (vector_length = block_size \* iidReps)

Example:

1.  Start from a vector CGF (e.g. independent Poissons),
2.  Apply a linear map `Y = A X`,
3.  Replicate the resulting vector-valued observation `Y` across blocks.

This makes `vector_length` larger than 1:
`vector_length = block_size * iidReps`.

Below, `X` is a 5-vector of independent Poissons, `A` is `3 × 5`, so one
observation `Y` lives in `R^3`. We then create `B = 9` i.i.d. copies of
`Y`, so `vector_length = 3 * 9 = 27`.

``` r


set.seed(1)

lambda5 <- c(1, 2, 3, 4, 5)

A <- matrix(
  c(
    1, 0, 0, 0, 0,
    0, 1, 1, 0, 0,
    0, 0, 0, 1, 1
  ),
  nrow = 3,
  byrow = TRUE
)

B <- 9

# B replicates of a 3-vector observation: Y = A X
# vector_length will be B*3
cg_Y <- linearlyMappedCGF(PoissonCGF, A, iidReps = B)

vector_length <- nrow(A) * B
n_sims <- 20000

Y <- cg_Y$rsim(n = n_sims, vector_length = vector_length, parameter_vector = lambda5)
dim(Y)
#> [1]    27 20000
```

A simple numerical check (untilted, `tvec = 0`): the mean of `Y = A X`
should be `E[Y] = A E[X] = A lambda`.

``` r


mu_Y <- as.vector(A %*% lambda5)

# Row means for the 27-vector, reshaped into 3 x B block layout:
mu_hat_blocks <- matrix(rowMeans(Y), nrow = nrow(A), ncol = B)

# Each column (each block) should be close to mu_Y:
mu_hat_blocks
#>         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]   [,7]    [,8]   [,9]
#> [1,] 0.99270 0.99730 1.00150 1.00280 1.01035 1.00300 0.9799 0.99150 0.9976
#> [2,] 4.98965 5.01495 4.98525 5.02030 5.01635 4.99945 5.0026 5.01825 4.9985
#> [3,] 9.03350 9.00030 8.96525 8.96055 9.01445 8.98030 9.0451 8.99960 8.9515
mu_Y
#> [1] 1 5 9
```

The same idea under tilting: If you tilt `Y` by a block tilt `h` (length
3), then `X` is implicitly tilted by `A^T h`, so under the tilted law:

- `E[X] = lambda * exp(A^T h)` (elementwise),
- `E[Y] = A E[X]`.

``` r


h_block <- c(0.2, -0.1, 0.15)
h_full  <- rep(h_block, times = B)

Y_tilt <- cg_Y$rsim(
  n = n_sims,
  vector_length = vector_length,
  parameter_vector = lambda5,
  tvec = h_full
)

lambda_tilt <- lambda5 * exp(as.vector(t(A) %*% h_block))
mu_Y_tilt <- A %*% lambda_tilt

mu_hat_tilt_blocks <- matrix(rowMeans(Y_tilt), nrow = nrow(A), ncol = B)

mu_hat_tilt_blocks
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,]  1.23320  1.21580  1.22385  1.21360  1.21035  1.21320  1.22380  1.22580
#> [2,]  4.50050  4.53055  4.50610  4.51885  4.52400  4.50105  4.53535  4.52110
#> [3,] 10.45365 10.46560 10.46960 10.44355 10.41785 10.46100 10.46205 10.45945
#>          [,9]
#> [1,]  1.22195
#> [2,]  4.51240
#> [3,] 10.50820
mu_Y_tilt
#>           [,1]
#> [1,]  1.221403
#> [2,]  4.524187
#> [3,] 10.456508
```

A quick “large vector” example:

``` r


set.seed(1)

d <- 2000
mu    <- rep(0, d)
sigma <- rep(1, d)
param <- c(mu, sigma)  # NormalCGF uses c(mu[1:d], sigma[1:d])

Xn <- NormalCGF$rsim(n = 2, vector_length = d, parameter_vector = param)
dim(Xn)
#> [1] 2000    2
```

## Exponential tilting via `tvec`

We now support tilted simulation by passing a nonzero `tvec` to
`rsim()`.

Informally, simulating with `tvec = h` means simulating from the
exponentially tilted law with density (or pmf) proportional to
`exp(h^T x)`.

Equivalently, if `K(t)` is the original CGF, the tilted CGF is

``` math
K_h(t) = K(t + h) - K(h).
```

For common exponential-family models, the tilted law stays in the same
family, often with a simple parameter update.

Not every `tvec` is valid: when the CGF only exists on a domain
(e.g. `t < rate` for Exponential/Gamma), tilted simulation requires
`tvec` to lie in that domain.

### Poisson

$`\lambda_h = \lambda e^{h}`$

``` r


set.seed(1)

B <- 20000
lambda <- 3
h <- 0.4

Xh <- PoissonCGF$rsim(
  n = B,
  vector_length = 1,
  parameter_vector = lambda,
  tvec = h
)

mean(Xh)
#> [1] 4.4778
lambda * exp(h)
#> [1] 4.475474
```

### Multinomial

Here the tilt is vector-valued: `tvec` has length `d`.

``` r


set.seed(1)

B <- 3000
d <- 30
N <- 50

# Random probability vector
p <- stats::runif(d)
p <- p / sum(p)

# A small tilt
h <- stats::rnorm(d, sd = 0.25)

p_h <- p * exp(h)
p_h <- p_h / sum(p_h)

Xh_multi <- MultinomialCGF$rsim(
  n = B,
  vector_length = d,
  parameter_vector = c(N, p),
  tvec = h
)

# Average proportions across simulations should be close to p_h
prop_hat <- rowMeans(Xh_multi) / N
c(
  rms_error = sqrt(mean((prop_hat - p_h)^2)),
  max_abs_error = max(abs(prop_hat - p_h))
)
#>     rms_error max_abs_error 
#>  0.0004134259  0.0007863824

head(cbind(prop_hat = prop_hat, p_h = p_h), 8)
#>        prop_hat        p_h
#> [1,] 0.01631333 0.01654760
#> [2,] 0.02311333 0.02335956
#> [3,] 0.04641333 0.04571452
#> [4,] 0.07073333 0.07028831
#> [5,] 0.01510667 0.01474633
#> [6,] 0.07135333 0.07124861
#> [7,] 0.07161333 0.07239972
#> [8,] 0.04263333 0.04243257
```

## A compositional example: sum of independent components

[`sumOfIndependentCGF()`](https://godrick.github.io/saddlepoint/reference/sumOfIndependentCGF.md)
constructs a CGF for a sum of independent random vectors. Simulation
composes in the obvious way: simulate each component and add.

``` r


set.seed(1)

# Y = Y1 + Y2, with Y1 ~ Poisson(2), Y2 ~ Poisson(5)
cg1 <- PoissonModelCGF(lambda = function(theta) theta[1])
cg2 <- PoissonModelCGF(lambda = function(theta) theta[2])

cg_sum <- sumOfIndependentCGF(list(cg1, cg2))

theta <- c(2, 5)
Y <- cg_sum$rsim(n = 10000, vector_length = 1, parameter_vector = theta)

mean(Y) # should be close to 7
#> [1] 7.0052
```
