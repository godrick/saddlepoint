# Creating custom CGFs (from vectorized functions)

This vignette shows how to build a new `CGF` object when you can write
down:

- the CGF $`K(t;\theta)`$ and its derivatives (at least up to order 4),
  and
- (optionally) the domain constraint(s) for $`t`$.

The recommended constructor for most univariate or elementwise models
is:

- [`createCGFfromVectorizedFunctions()`](https://godrick.github.io/saddlepoint/reference/createCGFfromVectorizedFunctions.md)

where you provide vectorized derivative functions that operate
elementwise on `tvec`.

## What you must provide

For
[`createCGFfromVectorizedFunctions()`](https://godrick.github.io/saddlepoint/reference/createCGFfromVectorizedFunctions.md),
you provide 5 vectorized functions:

- `K_vectorized(tvec, param)` returning a numeric vector of length
  `length(tvec)`
- `K1_vectorized(tvec, param)` returning a numeric vector
- `K2_vectorized(tvec, param)` returning a numeric vector (diagonal of
  the Hessian)
- `K3_vectorized(tvec, param)` returning a numeric vector
- `K4_vectorized(tvec, param)` returning a numeric vector

Often very useful:

- `ineq_constraint(tvec, param)` returning a numeric vector of
  constraints $`g(t,\theta)`$ with feasibility defined by $`g \le 0`$,
- `analytic_tvec_hat(x,param)` if you can solve $`K_1(\hat t;\theta)=x`$
  in closed form

## Example:

Let $`X \sim \mathrm{N}(\mu,\sigma^2)`$. The CGF is
``` math
K(t;\mu,\sigma^2) = \mu t + \frac{1}{2}\sigma^2 t^2,
```
and the derivatives are
``` math
K_1(t) = \mu + \sigma^2 t,\quad
K_2(t) = \sigma^2,\quad
K_3(t)=0,\quad
K_4(t)=0.
```

Below we implement a Normal CGF parameterized by
$`(\mu, \log\sigma^2)`$.

``` r


NormalLogSigma2CGF <- createCGFfromVectorizedFunctions(
  K_vectorized = function(tvec, param) {
    mu <- param[1]
    s2 <- exp(param[2])
    mu * tvec + 0.5 * s2 * tvec * tvec
  },
  K1_vectorized = function(tvec, param) {
    mu <- param[1]
    s2 <- exp(param[2])
    mu + s2 * tvec
  },
  K2_vectorized = function(tvec, param) {
    s2 <- exp(param[2])
    rep(s2, length(tvec))
  },
  K3_vectorized = function(tvec, param) rep(0, length(tvec)) * param[1],
  K4_vectorized = function(tvec, param) rep(0, length(tvec)) * param[1],

  # Optional: analytic saddlepoint t-hat for K1(t)=x
  analytic_tvec_hat = function(x, param) {
    mu <- param[1]
    s2 <- exp(param[2])
    (x - mu) / s2
  },
  rsim = function(n, vector_length, parameter_vector, tvec = NULL, ...) {
    mu <- parameter_vector[1]
    s2 <- exp(parameter_vector[2])

    if (!is.finite(mu) || !is.finite(s2) || s2 <= 0) stop("NormalLogSigma2CGF$rsim: parameters must satisfy sigma2 > 0 and be finite.")
    
    if (is.null(tvec)) tvec <- rep(0, vector_length)

    # Exponential tilting for Normal:
    # If X ~ N(mu, sigma^2), then under tilt t the law is N(mu + sigma^2 t, sigma^2).
    mu_tilt <- mu + s2 * tvec
    sd <- sqrt(s2)

    matrix(
      stats::rnorm(
        n = n * vector_length,
        mean = rep.int(mu_tilt, times = n),
        sd   = sd
      ),
      nrow = vector_length,
      ncol = n
    )
  },


  op_name = "NormalLogSigma2CGF"
)

# Simulate iid data (using the CGF's rsim method)
set.seed(1)
B <- 50
mu_true <- 1.5
s2_true <- 2.0
theta_true <- c(mu_true, log(s2_true))

# iid CGF for (Y_1,...,Y_B)
cgB <- iidReplicatesCGF(NormalLogSigma2CGF, iidReps = B, block_size = 1)

# Simulate one dataset (length B)
y <- as.numeric(cgB$rsim(
  n = 1,
  vector_length = B,
  parameter_vector = theta_true,
  flatten = TRUE
))
# Saddlepoint MLE (should match exact normal MLE very closely; for Normal SPA is exact)
fit <- find.saddlepoint.MLE(
  observed.data  = y,
  cgf            = cgB,
  starting.theta = c(0.01, log(1.3)),
  method         = "two_step"
)

# Exact MLE under Normal(mu, sigma^2):
mu_hat_exact <- mean(y)
s2_hat_exact <- mean((y - mu_hat_exact)^2)

cbind(
  spa_mle   = fit$MLEs.theta,
  exact_mle = c(mu_hat_exact, log(s2_hat_exact)),
  true      = c(mu_true, log(s2_true))
)
#>        spa_mle exact_mle      true
#> [1,] 1.6420553 1.6420553 1.5000000
#> [2,] 0.3036414 0.3036414 0.6931472
```

## Adding a domain constraint (when needed)

Some distributions have a restricted CGF domain in $`t`$. If
$`K(t;\theta)`$ is only defined when $`g(t,\theta)\le 0`$, implement:

``` r


ineq_constraint = function(tvec, param) {
  # return a numeric vector; feasible if all entries <= 0
}
```

Example: if a univariate CGF is defined for $`t < b(\theta)`$, you can
return `tvec - b(theta)`.
