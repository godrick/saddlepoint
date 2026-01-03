# CGF for a randomly-stopped sum

Builds a CGF object for the random vector \$\$Y = \sum\_{i=1}^{N}
X_i,\$\$ where:

- \\N\\ is a non-negative integer-valued scalar random variable with CGF
  `count_cgf`,

- \\X_i\\ are i.i.d. copies of a (possibly vector-valued) summand with
  CGF `summand_cgf`,

- \\N\\ is independent of all \\X_i\\.

## Usage

``` r
randomlyStoppedSumCGF(
  count_cgf,
  summand_cgf,
  block_size = NULL,
  iidReps = NULL,
  ...
)
```

## Arguments

- count_cgf:

  A `CGF` object for the (scalar) count variable \\N\\.

- summand_cgf:

  A `CGF` object for the summand \\X\\.

- block_size:

  Optional. The dimension \\d\\ of one observation of \\Y\\ (and hence
  of \\X\\). This can be a positive integer (fixed \\d\\) or `NULL`

- iidReps:

  Optional. Replication count \\B\\ for i.i.d. observations
  \\Y_1,\ldots,Y_B\\. May be a positive integer, `"any"`, or `NULL`.

- ...:

  Additional named arguments passed to CGF creation.

## Value

A `CGF` object.

## Details

**What do `block_size` and `iidReps` mean for RSS models?**

The replication arguments refer to i.i.d. replication of the variable
\\Y\\, not replication of summands inside the random sum.

If you have \\B\\ independent observations, \$\$Y_b = \sum\_{i=1}^{N_b}
X\_{b,i}, \qquad b=1,\ldots,B,\$\$ then:

- `block_size` is the dimension of one \\Y_b\\ (which equals the
  dimension of one \\X\_{b,i}\\).

- `iidReps` is the number of independent observations \\B\\.

**Why do we require `block_size` and/or `iidReps`?**

Many summand CGFs in this package are vectorized: if you accidentally
feed a long vector `tvec` into an 'untagged' RSS CGF, the summand CGF
may interpret that as a higher-dimensional summand \\X\\ (or multiple
i.i.d. summands inside one \\X\\), which silently changes the model.

To prevent this, `randomlyStoppedSumCGF()` requires at least one of
`block_size` or `iidReps`.

The semantics are:

- If `iidReps` is `NULL` (default), then `block_size` must be provided
  and we set `iidReps="any"`. This means the number of observations is
  inferred at evaluation time as `B = length(tvec)/block_size`.

- If `iidReps` is `"any"`, then `block_size` must be provided.

- If `iidReps` is a positive integer, then `block_size` may be `NULL`
  (it will be inferred as `length(tvec)/iidReps`), but providing
  `block_size` explicitly is encouraged.

## Examples

``` r
if (FALSE) { # \dontrun{
## Scalar RSS, B observations:
K.N <- GeometricModelCGF(prob = adaptor(fixed_param = 0.4))
K.X <- BinomialModelCGF(n = adaptor(fixed_param = 1), p = adaptor(indices = 1))

# Either specify B explicitly...
K.U1 <- randomlyStoppedSumCGF(K.N, K.X, iidReps = 40)      # block_size inferred as length(tvec)/40

# ...or specify the scalar block size and let B be inferred from length(tvec):
K.U2 <- randomlyStoppedSumCGF(K.N, K.X, block_size = 1)    # iidReps defaults to "any"

## Vector RSS: summand is d-dimensional
d <- 3
K.Xv <- MultinomialModelCGF(n = adaptor(fixed_param = 1),
                            prob_vec = adaptor(indices = 1:d))
K.Y  <- randomlyStoppedSumCGF(PoissonCGF, K.Xv, block_size = d, iidReps = "any")
} # }
```
