# Concatenation of CGFs (concatenationCGF)

[`concatenationCGF()`](https://godrick.github.io/saddlepoint/reference/concatenationCGF.md)
builds a CGF for a concatenated random vector with independent
components.

If you have independent random vectors
``` math
X^{(1)} \in \mathbb{R}^{d_1},\;\dots,\; X^{(L)} \in \mathbb{R}^{d_L},
```
and define the concatenation
``` math
Y = \bigl(X^{(1)}, \dots, X^{(L)}\bigr) \in \mathbb{R}^{d_1+\cdots+d_L},
```
then the CGF factorizes additively over disjoint sub-blocks of $`t`$:
``` math
K_Y(t;\theta) = \sum_{i=1}^L K_i\!\left(t^{(i)};\theta\right),
\qquad
t = \bigl(t^{(1)},\dots,t^{(L)}\bigr).
```

This differs from
[`sumOfIndependentCGF()`](https://godrick.github.io/saddlepoint/reference/sumOfIndependentCGF.md),
where all summands share the same `tvec` and are added as
scalars/vectors.

## Replication semantics

For concatenation, the natural “one observation” dimension is fixed:

- `block_size = sum(component_dims)`.
- You only need to specify `iidReps` (either a positive integer or
  `"any"`).

## Example 1: a bivariate Poisson model (exact likelihood available)

We observe $`B`$ i.i.d. bivariate vectors:
``` math
Y_b = (Y_{1b}, Y_{2b}), \qquad
Y_{1b} \sim \text{Poisson}(\lambda_1),\;
Y_{2b} \sim \text{Poisson}(\lambda_2),
```
independent across coordinates and across $`b`$.

The exact MLE is simply:
``` math
\hat\lambda_1 = \overline{Y}_1,\qquad \hat\lambda_2 = \overline{Y}_2.
```

``` r


set.seed(123)
B <- 30
lambda_true <- c(7, 15)

# Store data as a 2 x B matrix: columns are iid replicates
Y <- rbind(rpois(B, lambda_true[1]),
           rpois(B, lambda_true[2]))

# Model CGFs:
# lambda(theta) = theta[i], i=1,2  (so theta must be positive)
cg_1 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
cg_2 <- PoissonModelCGF(lambda = adaptor(indices = 2), iidReps = 1)

# Concatenate into a 2D vector, then replicate across columns (iidReps="any"):
cg_pp <- concatenationCGF(
  list(cg_1, cg_2),
  component_dims = 1,
  iidReps = "any"
)

# Saddlepoint MLE
fit <- find.saddlepoint.MLE(
  observed.data  = Y,                
  cgf            = cg_pp,
  starting.theta = c(1, 1.2),
  lb.theta       = c(1e-8, 1e-3),
  method         = "two_step"
)


cbind(
  true   = lambda_true,
  exact_mle  = rowMeans(Y),
  spa_mle = fit$MLEs.theta
)
#>      true exact_mle   spa_mle
#> [1,]    7   7.60000  7.600001
#> [2,]   15  13.86667 13.866665
```

Use
[`concatenationCGF()`](https://godrick.github.io/saddlepoint/reference/concatenationCGF.md)
when you want to model a multivariate observation by stitching together
independent subcomponents whose CGF objects already exist.
