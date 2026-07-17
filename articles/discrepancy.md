# Discrepancy diagnostic for saddlepoint MLEs

[`find.saddlepoint.MLE()`](https://godrick.github.io/saddlepoint/reference/find.saddlepoint.MLE.md)
maximizes a saddlepoint approximation to the log-likelihood. When
`discrepancy = TRUE`, it also returns a discrepancy approximation
intended to quantify how far the saddlepoint MLE is from the (possibly
unknown) exact MLE.

## What `discrepancy` means

Let

- $`\hat\theta_{\text{spa}}`$ be the saddlepoint MLE (from maximizing
  the saddlepoint likelihood),
- $`\hat\theta_{\text{true}}`$ be the exact MLE (from maximizing the
  exact likelihood).

The **true discrepancy** is
``` math
\delta = \hat\theta_{\text{true}} - \hat\theta_{\text{spa}}.
```

In practice $`\hat\theta_{\text{true}}`$ is unknown, so the package
returns an **approximated discrepancy** $`\hat\delta`$. A alternative
way to view it is:

- `theta_adj = theta_spa + discrepancy` is a second-order adjusted
  estimate.
- Comparing `abs(discrepancy) / std.error` gives a quick “is this
  difference practically large?” check.

## How to request it

``` r


fit <- find.saddlepoint.MLE(
  observed.data = y,
  cgf           = cgf,
  starting.theta = theta0,
  lb.theta       = ...,
  ub.theta       = ...,
  method         = "two_step",
  std.error      = TRUE,
  discrepancy    = TRUE
)

fit$MLEs.theta
fit$discrepancy
theta_adj <- fit$MLEs.theta + fit$discrepancy
```

## Example: Gamma(shape = alpha, rate = 1) with unknown shape (1 parameter)

``` r


fixed_rate <- 1
cgf_alpha <- GammaModelCGF(shape = adaptor(indices = 1), rate = adaptor(fixed_param = fixed_rate))


# simulate data 
set.seed(1)
n <- 80
alpha_true <- 6
x <- as.numeric(cgf_alpha$rsim(n = n, vector_length = 1, parameter_vector = alpha_true))

# true MLE
alpha_mle_true_fn <- function(alpha_, dat) -sum(dgamma(x = dat, shape = alpha_, 
                                                       rate = fixed_rate, log = TRUE))
alpha_true_mle <- nlminb(start = 0.3, objective = alpha_mle_true_fn, lower = 1e-03, dat = x)$par

# saddlepoint MLE + discrepancy 
fit <- find.saddlepoint.MLE(
  observed.data   = x,
  cgf             = cgf_alpha,
  starting.theta  = 2,
  lb.theta        = 1e-8,
  method          = "two_step",
  std.error       = TRUE,
  discrepancy     = TRUE
)

alpha_spa_mle <- fit$MLEs.theta
delta_hat <- fit$discrepancy

c(alpha_true = alpha_true,
  alpha_true_mle = alpha_true_mle,
  alpha_spa_mle   = alpha_spa_mle)
#>     alpha_true alpha_true_mle  alpha_spa_mle 
#>       6.000000       6.040570       6.027878

alpha_true_mle - alpha_spa_mle
#> [1] 0.01269219
delta_hat
#> [1] 0.01276576
```

### True vs approximated discrepancy

In practice you do not know the true discrepancy $`\delta`$ (because you
do not know $`\hat\theta_{\text{true}}`$). Here we compute it only to
highlight the accuracy of the approximation in this controlled toy
setting.

``` r


set.seed(2)

B <- 40 # number of simulated datasets
n <- 80
alpha_true <- 6

alpha_true_mle_vec <- numeric(B)
alpha_spa_mle_vec  <- numeric(B)
delta_true_vec  <- numeric(B)
delta_hat_vec   <- numeric(B)

for (b in 1:B) {
  x <- as.numeric(cgf_alpha$rsim(n = n, vector_length = 1, parameter_vector = alpha_true))

  # true MLE
  alpha_true_mle <- nlminb(start = 0.3, objective = alpha_mle_true_fn, lower = 1e-03, dat = x)$par

  # saddlepoint MLE + discrepancy
  fit <- find.saddlepoint.MLE(
    observed.data   = x,
    cgf = cgf_alpha,
    starting.theta  = 2,
    lb.theta = 1e-8,
    method  = "two_step",
    discrepancy     = TRUE
  )

  alpha_true_mle_vec[b] <- alpha_true_mle
  alpha_spa_mle_vec[b]   <- fit$MLEs.theta
  delta_true_vec[b]  <- alpha_true_mle - fit$MLEs.theta
  delta_hat_vec[b]   <- fit$discrepancy
}

op <- par(mfrow = c(1, 2))

plot(alpha_spa_mle_vec, alpha_true_mle_vec,
     xlab = "spa MLE ",
     ylab = "true MLE ")
abline(0, 1)

plot(delta_hat_vec, delta_true_vec,
     xlab = "Approximated discrepancy",
     ylab = "True discrepancy")
abline(0, 1)
```

![Two scatterplots comparing SPA vs exact MLE and estimated vs true
discrepancy (diagonal reference
lines).](discrepancy_files/figure-html/unnamed-chunk-3-1.png)

``` r


par(op)
```

### Notes

- If the model is an exponential family in the relevant sense, the true
  discrepancy is known to be $`0`$. In such cases the discrepancy
  approximation will also return $`0`$.
- `discrepancy` can be useful as a diagnostic and as a small correction
  (e.g. if it is tiny relative to standard errors).
