# Create a Parametric Subunitary Multinomial CGF Object

Creates a CGF object for the subunitary multinomial distribution with
total count N = n(theta) and active-category probabilities pi =
prob_vec(theta), where sum(pi) \<= 1.

## Usage

``` r
SubunitaryMultinomialModelCGF(n, prob_vec, iidReps = "any", ...)
```

## Arguments

- n:

  A function (or adaptor) mapping theta to the multinomial total count
  N.

- prob_vec:

  A function (or adaptor) mapping theta to the active-category
  probability vector pi (length d). These entries should be non-negative
  and sum(pi) should be \<= 1.

- iidReps:

  Either "any" or a positive integer specifying how many i.i.d. blocks
  are expected.

- ...:

  Additional named arguments passed to CGF creating functions.

## Value

A CGF object.

## Details

Let \\Y \sim \mathrm{Multinomial}(N,\pi_1,\dots,\pi_d,\pi\_{d+1})\\ with
\\\sum\_{i=1}^{d+1} \pi_i = 1\\. Conditioning on \\\\Y\_{d+1}=0\\\\, we
have \\W = (Y \mid Y\_{d+1}=0) \sim
\mathrm{Multinomial}(N,p_1,\dots,p_d) \\, where \$\$ p_i =
\frac{\pi_i}{\sum\_{j=1}^d \pi_j}, \quad i=1,\dots,d. \$\$

**The CGF:**

From the perspective of restricting to \\\\Y\_{d+1}=0\\\\, one
effectively sets \\t\_{d+1}=-\infty\\ in the multinomial CGF. This
yields the *subunitary* CGF: \$\$
K\_{\mathrm{subunitary}}(t_1,\dots,t_d) =
\log\bigl\\\Pr(Y\_{d+1}=0)\bigr\\ \\+\\ K\_{W}(t_1,\dots,t_d)\\, \$\$
where \\\log\bigl\\\Pr(Y\_{d+1}=0)\bigr\\ = N \log\left(\\\sum\_{i=1}^d
\pi_i\right)\\ and \\K\_{W}\\ is the usual multinomial CGF for
\\\mathrm{Multinomial}(N,p_1,\dots,p_d)\\.

In many applications (such as certain capture-recapture models), some
categories are effectively impossible based on observed data, so they
can be merged into a single "impossible" category with zero count. The
probabilities of the \\d\\ active categories then sum to less than one.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
d <- 3
N <- 10
pi_vec <- c(0.2, 0.3, 0.1)  # sum < 1
theta <- c(N, pi_vec)

# Build a theta-mapped CGF: theta -> (N, pi_1,...,pi_d)
cgf <- SubunitaryMultinomialModelCGF(
  n        = function(th) th[1],
  prob_vec = function(th) th[2:(d + 1)],
  iidReps  = 2
)

# With iidReps = 2, tvec must be length 2*d:
t_block <- rnorm(d)
tvec <- c(t_block, t_block)

# Evaluate K and K1 at theta:
cgf$K(tvec, theta)
cgf$K1(tvec, theta)
} # }
```
