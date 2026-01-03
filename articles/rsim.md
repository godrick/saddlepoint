# Random simulation with CGF objects: rsim

This vignette documents how `rsim(...)` works across CGF objects and
wrappers in this codebase. We will keep extending it as new CGFs gain
simulation support.

> **Conventions used here**
>
> - In cgf\$rsim(iidReps = B, …), iidReps = number of simulated draws.
> - In building blocks like iidReplicatesCGF(…, iidReps = …), iidReps =
>   how many i.i.d. blocks tvec represents (or “any” to infer).
> - For a `d`-dimensional random vector, simulation returns a
>   `d x (iidReps)` matrix where each **column** is one draw.
> - If `d == 1` and `drop = TRUE` (default), simulation returns a
>   numeric vector of length `iidReps`.

## What `block_size` is, and when you need to think about it

A lot of CGFs in this `saddlepoint` package are *vectorized*. That means
a long `tvec` could mean either:

- **one** high-dimensional random vector, or
- **B** independent “blocks” of a smaller random vector (i.i.d.
  replication).

The wrapper
[`iidReplicatesCGF()`](https://godrick.github.io/saddlepoint/reference/iidReplicatesCGF.md)
disambiguates those interpretations. It uses:

- `block_size` = dimension of one observation (one block), and
- `iidReps` (inside the wrapper) = how many blocks to expect (or `"any"`
  to infer it from `length(tvec)`).

### Does `block_size` matter for `rsim()`?

Usually, **no**. `rsim(iidReps = B, ...)` is unambiguous: you asked for
`B` draws, and simulation can determine the dimension `d` from the model
parameters.

You mainly need `block_size` when *constructing* a CGF (e.g., RSS
models) so that evaluation methods (`K`, `K1`, `K2`, etc.) cannot
silently reinterpret a long `tvec` as “a different statistical model”.

Wrappers like
[`linearlyMappedCGF()`](https://godrick.github.io/saddlepoint/reference/linearlyMappedCGF.md)
can usually infer `block_size` internally (e.g., `nrow(A)`), so you
rarely worry about it there.

## Quick checklist for `rsim()` support

A CGF supports simulation if:

``` r
# cgf$has_simulate()
```

returns `TRUE`.

When it does, the canonical call is:

``` r
Y <- cgf$rsim(iidReps = 10, parameter_vector = theta)
```

## Example 1a: PoissonCGF (scalar and vector)

``` r
set.seed(1)

# Scalar Poisson (d = 1)
lambda <- 2
x <- PoissonCGF$rsim(iidReps = 10, parameter_vector = lambda)
x
#>  [1] 1 1 2 4 1 4 4 2 2 0

# Vector Poisson (d = 2) interpreted as independent components
lambda2 <- c(2, 3)
X <- PoissonCGF$rsim(iidReps = 5, parameter_vector = lambda2, drop = FALSE)
X
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    3    3    3    1
#> [2,]    1    2    3    8    4
```

## Example 1b: NormalCGF (vectorized independent normals)

``` r
set.seed(1)

# Two independent normals packed into one CGF parameter vector:
# parameter_vector = c(mu1, mu2, sigma1, sigma2)
param <- c(0, 5, 1, 2)

X <- NormalCGF$rsim(iidReps = 5, parameter_vector = param, drop = FALSE)
X
#>            [,1]       [,2]      [,3]      [,4]      [,5]
#> [1,] -0.6264538 -0.8356286 0.3295078 0.4874291 0.5757814
#> [2,]  5.3672866  8.1905616 3.3590632 6.4766494 4.3892232
dim(X)  # 2 x 5
#> [1] 2 5
rowMeans(X)  # ~ c(0, 5)
#> [1] -0.01387285  5.55655682
```

## Example 1c: BinomialCGF (vectorized)

``` r
set.seed(1)

# Two independent binomials:
# parameter_vector = c(n1, n2, p1, p2)
param <- c(10, 20, 0.2, 0.6)

X <- BinomialCGF$rsim(iidReps = 6, parameter_vector = param, drop = FALSE)
X
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    2    1    4    2    1
#> [2,]   13    9    9   11   15   14
rowMeans(X)              # ~ c(2, 12)
#> [1]  1.833333 11.833333
rowMeans(X) / c(10, 20)  # ~ c(0.2, 0.6)
#> [1] 0.1833333 0.5916667
```

## Example 2: linearlyMappedCGF simulation (sum of independent Poissons)

Let `X = (X1, X2)` with `X1 ~ Pois(2)`, `X2 ~ Pois(3)` independent, and
define `Y = X1 + X2`. This is implemented as a linear map with
`A = [1 1]`.

``` r
set.seed(1)

cg_vec <- PoissonCGF
A <- matrix(c(1, 1), nrow = 1)

cg_sum <- linearlyMappedCGF(cg_vec, A)

lambda <- c(2, 3)
Y <- cg_sum$rsim(iidReps = 10000, parameter_vector = lambda)

mean(Y)         # should be ~ 5
#> [1] 5.0059
var(Y)          # should be ~ 5
#> [1] 4.958761
```

## Example 3: Parameter mapping with adaptors (avoiding “wrong dimension” mistakes)

A common pitfall is passing a combined parameter vector to a component
CGF that expects only one number. For example, if `count_cgf` is
`PoissonCGF` and you pass `theta = c(lambdaN, lambdaX)`, then
`PoissonCGF` interprets this as a **2-dimensional** Poisson parameter,
not a scalar.

The fix is to *adapt* each sub-CGF so it receives only the piece of
`theta` it needs.

``` r
set.seed(1)

# theta = c(lambdaN, lambdaX)
theta <- c(lambdaN = 2, lambdaX = 3)

count_cgf   <- PoissonModelCGF(lambda = adaptor(indices = 1))
summand_cgf <- PoissonModelCGF(lambda = adaptor(indices = 2))

rss <- randomlyStoppedSumCGF(count_cgf, summand_cgf, block_size = 1)

Y <- rss$rsim(iidReps = 10000, parameter_vector = theta)
mean(Y)  # for compound Poisson, E[Y] = E[N] * E[X] = lambdaN * lambdaX = 6
#> [1] 5.9946
```

## Example 4: RSS with vector summands (Poisson thinning check)

Let `N ~ Poisson(lambda)` and `X_i ~ Multinomial(1, p)` a one-hot vector
(categorical draw). Then `Y = sum_{i=1}^N X_i` should have independent
Poisson components with means `lambda * p`.

``` r
set.seed(1)

d <- 3
theta <- c(lambda = 10, p = c(0.2, 0.3, 0.5))

count_cgf <- PoissonModelCGF(lambda = adaptor(indices = 1))
summand_cgf <- MultinomialModelCGF(
  n        = adaptor(fixed_param = 1),
  prob_vec = adaptor(indices = 2:(d + 1))
)

rss_vec <- randomlyStoppedSumCGF(count_cgf, summand_cgf, block_size = d)

Y <- rss_vec$rsim(iidReps = 20000, parameter_vector = theta, drop = FALSE)

rowMeans(Y)  # should be ~ lambda * p = (2, 3, 5)
#> [1] 2.00410 3.00845 5.01905
```

## Example 5: MultinomialCGF simulation (odds vs probabilities)

`MultinomialCGF` expects `parameter_vector = c(N, x_1, ..., x_d)`. The
`x` vector is treated as odds (or as probabilities if it already sums to
1).

``` r
set.seed(1)

N <- 10
odds <- c(2, 3, 5)
param_odds <- c(N, odds)

X <- MultinomialCGF$rsim(iidReps = 6, parameter_vector = param_odds, drop = FALSE)
X
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    2    1    4    2    1
#> [2,]    3    5    5    3    1    2
#> [3,]    6    3    4    3    7    7
colSums(X)  # should all equal N
#> [1] 10 10 10 10 10 10

# Same distribution if you pass probabilities instead of odds
p <- odds / sum(odds)
param_p <- c(N, p)

Xp <- MultinomialCGF$rsim(iidReps = 6, parameter_vector = param_p, drop = FALSE)
colSums(Xp)
#> [1] 10 10 10 10 10 10
```

## Example 6: MultinomialModelCGF (theta-mapped multinomial)

``` r
set.seed(1)

d <- 3
theta <- c(N = 12, p = c(0.1, 0.2, 0.7))

cgf <- MultinomialModelCGF(
  n        = function(th) th[1],
  prob_vec = function(th) th[2:(d + 1)],
  iidReps  = "any"
)

X <- cgf$rsim(iidReps = 5, parameter_vector = theta, drop = FALSE)
X
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    3    1
#> [2,]    2    4    5    2    0
#> [3,]   10    7    7    7   11
colSums(X)
#> [1] 12 12 12 12 12
```

## Example 7: SubunitaryMultinomialCGF

`SubunitaryMultinomialCGF` uses parameters `c(N, pi_1,...,pi_d)` where
`sum(pi) <= 1` and conceptually corresponds to a (d+1)-category
multinomial **conditioned** on the last category being zero.

Simulation returns the `d` active-category counts and still sums to `N`.

``` r
set.seed(1)

N <- 10
pi <- c(0.2, 0.3, 0.1)  # sum < 1 (subunitary)
param <- c(N, pi)

W <- SubunitaryMultinomialCGF$rsim(iidReps = 6, parameter_vector = param, drop = FALSE)
W
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    2    4    2    6    4    2
#> [2,]    6    3    4    3    6    7
#> [3,]    2    3    4    1    0    1
colSums(W)
#> [1] 10 10 10 10 10 10
```

## Example 8: sumOfiidCGF (sum of i.i.d. Poissons)

``` r
set.seed(1)

base <- PoissonModelCGF(lambda = adaptor(indices = 1))
sum10 <- sumOfiidCGF(base, n = 10, block_size = 1)

theta <- c(lambda = 2)
Y <- sum10$rsim(iidReps = 20000, parameter_vector = theta)

mean(Y)  # should be ~ 20
#> [1] 19.99555
var(Y)   # should be ~ 20
#> [1] 19.69872
```

## Example: sumOfIndependentCGF (sum of independent, not-necessarily-identical terms)

``` r
set.seed(1)

theta <- c(lambda1 = 2, lambda2 = 3)

cg1 <- PoissonModelCGF(lambda = adaptor(indices = 1))
cg2 <- PoissonModelCGF(lambda = adaptor(indices = 2))

cg_sum <- sumOfIndependentCGF(list(cg1, cg2), block_size = 1)

Y <- cg_sum$rsim(iidReps = 20000, parameter_vector = theta)

mean(Y)  # should be ~ 5
#> [1] 4.9943
var(Y)   # should be ~ 5
#> [1] 4.993817
```

## Other issues

### 1) RSS: `count_cgf` must be scalar

In RSS models, `count_cgf` represents a *scalar* count `N`, so
simulation requires a scalar draw per replicate. If your count CGF is
vectorized (like `PoissonCGF`) and you accidentally pass a vector
parameter, you’ll get a “count must be scalar” error.

The fix is to use parameter adaptors (as shown above) so the count CGF
receives only the scalar parameter it needs.

------------------------------------------------------------------------

## Adding a new `rsim()` implementation

For any new base distribution, define a simulator and pass it to
`createCGF(rsim = ...)`.

**Contract expected by `CGF$rsim`:**

- simulator returns either:
  - a numeric vector of length `d * iidReps`, or
  - a numeric `d x iidReps` matrix.
- each column corresponds to one simulated draw of the underlying random
  vector.

Then wrappers can often be implemented in a small, compositional way:
e.g. `linearlyMappedCGF` simulates `X` and returns `A %*% X`.
