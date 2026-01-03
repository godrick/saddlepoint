# CGF Object of a linearly mapped random variable \\Y = A \\ X\\

Creates a CGF object for the random vector \\Y = A(\theta) \\ X\\, where
\\X\\ is described by the input CGF `cgf`. The argument `matrix_A` can
be:

- **A numeric matrix** (dense or sparse).

- **A function**: \\\theta \mapsto A(\theta)\\ See details for the
  possible options of `matrix_A`.

If `matrix_A` is a function, it is called for each invocation of the CGF
methods to retrieve the current matrix (allowing parameter-dependent
transformations).

## Usage

``` r
linearlyMappedCGF(cgf, matrix_A, iidReps = "any", ...)
```

## Arguments

- cgf:

  An object of class `CGF` for the base distribution \\X\\.

- matrix_A:

  Either a numeric matrix (dense or `Matrix` sparse), or a function
  `function(theta) -> A(theta)` returning one of: numeric dense matrix,
  `Matrix` sparse matrix, RTMB AD‑dense (an `advector` with a `dim`) or
  RTMB `adsparse`.

- iidReps:

  Either `"any"` (default) or a positive integer. See
  [`iidReplicatesCGF`](https://godrick.github.io/saddlepoint/reference/iidReplicatesCGF.md)
  for the replication semantics.

- ...:

  Additional named arguments passed to `CGF` creation functions.

## Value

A `CGF` object for \\Y = A \\ X\\.

## Details

Accepted types for `matrix_A`:

- **numeric constant** matrix, dense or `Matrix` sparse. If dense, it is
  converted once at construction to a sparse `Matrix`.

- **function** \\\theta \mapsto A(\theta)\\ returning one of:

  - a numeric `matrix` (dense); it will be internally handled as sparse;

  - a `Matrix` sparse matrix (preferred for large problems);

  - an RTMB::AD-dense matrix (RTMB `advector` with a `dim` attribute);

  - an RTMB::AD-sparse matrix (RTMB `adsparse`), created e.g. via
    `A = RTMB::AD(Matrix::sparseMatrix(...)); A@x[] = ...`.

## Examples

``` r
## Example 1: constant numeric A
if (FALSE) { # \dontrun{
lambda_fun <- function(theta) c(theta[1], theta[2])   # two Poisson rates
pois2 <- PoissonModelCGF(lambda = lambda_fun, iidReps = "any")

A_const <- rbind(c(1, 0),
                 c(0.5, 1))
mapped  <- linearlyMappedCGF(cgf = pois2, matrix_A = A_const, iidReps = "any")

B <- 3L # 3 replicated 2-d blocks
t_one <- c(0.10, -0.05)
y_one <- c(3.0,   4.0)
tvec  <- rep(t_one, B)
y <- rep(y_one, B)
theta <- c(2.0, 3.0)

## Identity check: K1_Y(t) = A K1_X(A^T t)
t_block <- t_one
A_now   <- A_const
lhs_K1  <- mapped$K1(t_block, theta)
poisk1 <- pois2$K1(as.vector(t(A_now) %*% t_block), theta)
rhs_K1  <- A_now %*% poisk1

## SPA negative log-likelihood + gradient + Hessian
res <- compute.spa.negll(
  parameter_vector = theta,
  observed.data    = y,
  cgf              = mapped,
  gradient         = TRUE,
  hessian          = TRUE,
  tvec_source      = "newton",
  spa_method       = "standard"
)
} # }
## Example 2: A(theta) dense-vs-sparse (adsparse)
if (FALSE) { # \dontrun{
library(Matrix)

lambda_fun <- function(theta) c(theta[1], theta[2])   # base 2-d Poisson
pois2 <- PoissonModelCGF(lambda = lambda_fun, iidReps = "any")

## Dense A(theta) depending on theta[1]
A_theta_dense <- function(theta) {
  matrix(c(1, 0,
           0.5*theta[1], 1), 2, 2, byrow = TRUE)
}

## AD-sparse A(theta) with fixed sparsity pattern, AD-filled values
A_theta_sparse <- function(theta) {
  A0 <- sparseMatrix(i = c(1L, 2L, 2L),
                     j = c(1L, 1L, 2L),
                     x = c(1, 0, 1),
                     dims = c(2L, 2L))
  A  <- RTMB::AD(A0)
  A@x[] = c(1, 0.5*theta[1], 1)
  A
}

mapped_dense  <- linearlyMappedCGF(cgf = pois2, matrix_A = A_theta_dense,  iidReps = "any")
mapped_sparse <- linearlyMappedCGF(cgf = pois2, matrix_A = A_theta_sparse, iidReps = "any")

## Data: B = 3 identical 2-vectors
B     <- 3L
y_one <- c(3.0, 4.0)
y     <- rep(y_one, B)
theta <- c(2.0, 3.0) # only theta[1] affects A(theta)

nll_dense <- compute.spa.negll(theta, y, mapped_dense,
                               gradient=TRUE, hessian=TRUE,
                               tvec_source="newton", spa_method="standard")
nll_sparse <- compute.spa.negll(theta, y, mapped_sparse,
                                gradient=TRUE, hessian=TRUE,
                                tvec_source="newton", spa_method="standard")
} # }
```
