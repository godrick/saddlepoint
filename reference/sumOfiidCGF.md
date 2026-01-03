# Create a CGF object for the sum of `n` i.i.d. random variables.

Given a base CGF object `cgf` describing a random variable \\X\\, this
function returns a new CGF object for the sum of `n` i.i.d. copies of
\\X\\, i.e. \\Y = X_1 + \cdots + X_n\\.

## Usage

``` r
sumOfiidCGF(cgf, n, block_size = NULL, iidReps = NULL, ...)
```

## Arguments

- cgf:

  A `CGF` object describing the base random variable \\X\\.

- n:

  Either a positive scalar, or a function/adaptor mapping \\\theta
  \mapsto n(\theta)\\.

- block_size:

  Optional. Positive integer giving the dimension of one observation of
  \\Y\\.

- iidReps:

  Optional. Either `NULL`, `"any"`, or a positive integer.

- ...:

  Additional arguments passed to CGF creation.

## Value

A `CGF` object for \\Y = \sum\_{i=1}^{n(\theta)} X_i\\.

## Details

\$\$K_Y(t;\theta) = n(\theta)\\K_X(t;\theta).\$\$

For every order \\r \ge 1\\, \$\$K_Y^{(r)}(t;\theta) =
n(\theta)\\K_X^{(r)}(t;\theta).\$\$

- `n` is the number of summands inside each \\Y\\.

- `iidReps`/`block_size` (if supplied) describe i.i.d. replication of
  \\Y\\ across observations \\Y_1,\ldots,Y_B\\.

**Replication (iidReps / block_size):**

- If both `iidReps` and `block_size` are `NULL`, no replication is
  applied.

- If `block_size` is provided but `iidReps` is `NULL`, we set
  `iidReps = "any"` and infer the number of blocks from
  `length(tvec) / block_size` at evaluation time.

- If `iidReps = "any"`, then `block_size` must be provided.

- If `iidReps` is a positive integer, `block_size` may be `NULL`, though
  providing it is safer.

## Examples

``` r
if (FALSE) { # \dontrun{
base_cgf <- NormalCGF
sum5 <- sumOfiidCGF(base_cgf, n=5)
theta <- c(1.2, 0.8)
tvec  <- c(0.1, -0.2, 0.0)

all.equal(sum5$K(tvec, theta), 5*base_cgf$K(tvec, theta))
all.equal(sum5$K1(tvec, theta), 5*base_cgf$K1(tvec, theta))
} # }
```
