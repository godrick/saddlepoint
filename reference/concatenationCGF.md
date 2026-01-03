# Concatenation of CGF objects

Constructs a `CGF` object for a concatenated random vector whose
components are independent and each specified by its own `CGF` object.

If `cgf_list = list(C1, ..., CL)` and `component_dims = (d1,...,dL)`,
then `tvec` is partitioned as \$\$tvec = (t^{(1)}, ..., t^{(L)})
\quad\text{with}\quad length(t^{(i)}) = d_i,\$\$ and the resulting CGF
satisfies \$\$K(tvec;\theta) = \sum\_{i=1}^L K_i(t^{(i)};\theta).\$\$

## Usage

``` r
concatenationCGF(cgf_list, component_dims = 1L, iidReps = "any", ...)
```

## Arguments

- cgf_list:

  A non-empty list of `CGF` objects.

- component_dims:

  An integer vector of length `length(cgf_list)` giving the dimension of
  each component. If a single integer is supplied, it is recycled.

- iidReps:

  Either `"any"` (default) or a positive integer. See Details.

- ...:

  Additional arguments passed to
  [`createCGF()`](https://godrick.github.io/saddlepoint/reference/createCGF.md)
  (advanced use).

## Value

A `CGF` object.

## See also

[`iidReplicatesCGF`](https://godrick.github.io/saddlepoint/reference/iidReplicatesCGF.md),
[`sumOfIndependentCGF`](https://godrick.github.io/saddlepoint/reference/sumOfIndependentCGF.md)

## Examples

``` r
## univariate + univariate (Poisson + Gamma)
# Use adaptors so both components share a single global parameter vector:
#   theta = c(lambda, shape, rate)
cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
cg_gam  <- GammaModelCGF(shape = adaptor(indices = 2),
                         rate  = adaptor(indices = 3),
                         iidReps = 1)

# component_dims can be a scalar (recycled) or a vector:
cg <- concatenationCGF(list(cg_pois, cg_gam),
                       component_dims = 1L,
                       iidReps = "any")

theta <- c(lambda = 2, shape = 5, rate = 3)

# Two iid blocks; each block is (Poisson, Gamma) so block_size = 2.
tvec <- c(0.05, 0.10,
         -0.02, 0.20)

# K is additive across blocks and across components:
K_manual <- cg_pois$K(tvec[1], theta) + cg_gam$K(tvec[2], theta) +
            cg_pois$K(tvec[3], theta) + cg_gam$K(tvec[4], theta)
stopifnot(all.equal(cg$K(tvec, theta), K_manual))

# K1 is blockwise concatenation:
K1_manual <- c(cg_pois$K1(tvec[1], theta), cg_gam$K1(tvec[2], theta),
               cg_pois$K1(tvec[3], theta), cg_gam$K1(tvec[4], theta))
stopifnot(all.equal(as.numeric(cg$K1(tvec, theta)),
                    as.numeric(K1_manual)))

# Analytic t-hat is available if all children have it:
if (isTRUE(cg$has_analytic_tvec_hat())) {
  x <- c(12, 5, 8, 7)  # observed data (same layout as tvec)
  t_hat <- cg$analytic_tvec_hat(x, theta)
  t_hat_manual <- c(cg_pois$analytic_tvec_hat(x[1], theta),
                    cg_gam$analytic_tvec_hat(x[2], theta),
                    cg_pois$analytic_tvec_hat(x[3], theta),
                    cg_gam$analytic_tvec_hat(x[4], theta))
  stopifnot(all.equal(as.numeric(t_hat), as.numeric(t_hat_manual)))
}

## Multivariate + univariate (Multinomial + Poisson)
# Multinomial dimension d = 3; concatenate with a Poisson scalar -> block_size = 4.
cg_mult <- MultinomialModelCGF(
  n        = adaptor(fixed_param = 10),
  prob_vec = adaptor(indices = 2:4),
  iidReps  = 1
)
cg_pois2 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)

cg2 <- concatenationCGF(list(cg_mult, cg_pois2),
                        component_dims = c(3L, 1L),
                        iidReps = "any")

theta2 <- c(lambda = 4, p1 = 0.2, p2 = 0.3, p3 = 0.5)

# Two iid blocks, each of length 4:
tvec2 <- c(0.02, -0.01, 0.03,  0.10,
          -0.02,  0.00, 0.01, -0.05)

# Shape checks:
stopifnot(length(cg2$K1(tvec2, theta2)) == length(tvec2))
K2_full <- cg2$K2(tvec2, theta2)
stopifnot(all(dim(K2_full) == c(length(tvec2), length(tvec2))))


## find.saddlepoint.MLE on a concatenation model
# \donttest{
set.seed(123)
B <- 30
lambda_true <- c(5, 15)
Y <- rbind(rpois(B, lambda_true[1]),
           rpois(B, lambda_true[2]))

cg_1 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
cg_2 <- PoissonModelCGF(lambda = adaptor(indices = 2), iidReps = 1)
cg_pp <- concatenationCGF(list(cg_1, cg_2),
                          component_dims = 1L,
                          iidReps = "any")

fit <- find.saddlepoint.MLE(
  observed.data  = Y, # 2 x B matrix -> columns are iid blocks
  cgf            = cg_pp,
  starting.theta = c(0.1, 0.9),
  lb.theta       = c(0.001, 0.1),
  method         = "two_step"
)
#> Treating columns of 'observed.data' as i.i.d. replicate blocks.

stopifnot(max(abs(fit$MLEs.theta - rowMeans(Y))) < 1e-5)
# }
```
