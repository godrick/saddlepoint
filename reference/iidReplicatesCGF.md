# Replicate a CGF object over multiple i.i.d. blocks and/or with a fixed block size

Extends a given `CGF` object to handle multiple i.i.d. blocks. You can
specify:

- `iidReps` only,

- `block_size` only,

- or both `iidReps` and `block_size`.

## Usage

``` r
iidReplicatesCGF(cgf, iidReps = "any", block_size = NULL)
```

## Arguments

- cgf:

  A `CGF` object.

- iidReps:

  Either `"any"` (default) or a positive integer.

- block_size:

  Either NULL, a positive integer, or a function
  `function(param) -> positive integer`. When a function is provided,
  the block size is determined dynamically from the current parameter
  vector.

## Value

A `CGF` object

## Details

Let \\N = \mathrm{length(tvec)}\\, `d = block_size`, and `iidReps` the
number of blocks.

- If `iidReps = "any"` and `block_size = NULL`: pass-through (no
  splitting).

- If `iidReps = "any"` and `block_size = d`: require `N %% d == 0`, set
  `iidReps = N/d`.

- If `iidReps = m` and `block_size = NULL`: require `N %% m == 0`, set
  `d = N/m`, `iidReps = m`.

- If `iidReps = m` and `block_size = d`: require `N == d * m`, set
  `iidReps = m`.

## Examples

``` r
## Base CGF: univariate Poisson with lambda(theta) = theta[1]
pois <- PoissonModelCGF(lambda = function(th) th[1], iidReps = "any")
theta <- c(2)  # rate

## iidReps = "any", block_size = NULL
pass <- iidReplicatesCGF(pois, iidReps = "any", block_size = NULL)
identical(pass, pois)  # TRUE
#> [1] TRUE
pass$K1(0, theta)  # same as pois$K1(0, theta)
#> [1] 2

## block_size only
tvec <- c(0.1, -0.2, 0.0, 0.3)    # N = 4
agg_b <- iidReplicatesCGF(pois, block_size = 2)  # 2 blocks each of length 2
k1_b  <- agg_b$K1(tvec, theta)
k1_ref <- c(pois$K1(tvec[1:2], theta), pois$K1(tvec[3:4], theta))
all.equal(k1_b, k1_ref) # TRUE
#> [1] TRUE

## iidReps only: split tvec into exactly iidReps blocks (block size inferred)
t6 <- c(0.1, -0.2, 0.0, 0.3, 0.2, -0.1)  # N = 6
agg_m <- iidReplicatesCGF(pois, iidReps = 3)     # B = 3, d = N/B = 2
K2_m  <- agg_m$K2(t6, theta)
# Off-block covariances are zero (block diagonal)
is_zero <- function(M) max(abs(M)) < 1e-12
is_zero(K2_m[1:2, 3:4]) && is_zero(K2_m[1:2, 5:6]) && is_zero(K2_m[3:4, 5:6])
#> [1] TRUE

# Both iidReps and block_size: require N == d * B
agg_fix <- iidReplicatesCGF(pois, iidReps = 3, block_size = 2)  # N must be 6 here
all.equal(agg_fix$K1(t6, theta), agg_m$K1(t6, theta))  # same objects agg_fix/agg_m
#> [1] TRUE
if (FALSE) { # \dontrun{
  # Mismatch example (errors at evaluation, not at construction):
  bad <- iidReplicatesCGF(pois, iidReps = 3, block_size = 4)
  bad$K1(t6, theta)  # length(t6) = 6 != 3 * 4  ==> error
} # }


## Here the inner Poisson builder requires EXACTLY 2 inner replicates.
inner2 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 2)
t6 <- c(0.1, -0.2, 0.0, 0.3, 0.2, -0.1)  # N = 6
# With block_size = 2, each outer block (length d = 2) satisfies the inner rule.
outer <- iidReplicatesCGF(inner2, block_size = 2)  # B = 3 blocks; d = 2 per block
outer$K1(t6, theta)  # valid; each block passes inner iidReps = 2 check
#> [1] 2.210342 1.637462 2.000000 2.699718 2.442806 1.809675

## Non-identical setup
## lambda(theta) returns a vector of length 2; inner builder set to "any".
nonid <- PoissonModelCGF(lambda = function(th) c(th[1], 3*th[1]), iidReps = "any")
t12 <- rep(c(0.05, -0.10, 0.00, 0.20), 3)  # N = 12
# With iidReps = 3 (outer), we have 3 blocks each of size 4.
# The inner builder sees d = 4 with 2 inner replicates per block.
agg_nonid <- iidReplicatesCGF(nonid, iidReps = 3)
agg_nonid$K1(t12, theta)
#>  [1] 2.102542 5.429025 2.000000 7.328417 2.102542 5.429025 2.000000 7.328417
#>  [9] 2.102542 5.429025 2.000000 7.328417
```
