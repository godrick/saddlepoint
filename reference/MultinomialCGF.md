# MultinomialCGF

A ready-to-use CGF object for handling multiple i.i.d. multinomial
random vectors. The number of multinomial random vectors is determined
by the length of `tvec` and the length of the multinomial probability
vector. For user-specified parameter mappings, refer to
[`MultinomialModelCGF`](https://godrick.github.io/saddlepoint/reference/MultinomialModelCGF.md).

## Usage

``` r
MultinomialCGF
```

## Format

An object of class `CGF` (an R6 class), with standard methods `K`, `K1`,
`K2`, `K3operator`, `K4operator`, etc.

## Details

**Parameter Vector Interpretation**: The `parameter_vector` used when
calling methods like `K1(tvec, parameter_vector)` should be of the form
\\(n, x_1, \dots, x_d)\\, where:

- \\n\\ is the total count (a positive integer),

- \\x = (x_1, \dots, x_d)\\ is a vector of non-negative entries, not all
  zero.

The vector \\x\\ can be interpreted in two ways:

- **Odds**: \\x_i\\ are odds values. Probabilities are derived as \\p_i
  = x_i / \sum\_{j=1}^{d} x_j\\.

- **Probabilities**: If \\x\\ sums to 1, it is treated directly as the
  probability vector \\p\\.

**IID Assumption**: All multinomial random vectors processed by
`MultinomialCGF` are assumed to be i.i.d.

## Examples

``` r
tvec <- rep(0,3)
parameter_vector <- c(10, 2, 3, 5)
# total count=10, odds=(2,3,5) or probabilities=c(0.2,0.3,0.5)
MultinomialCGF$K1(tvec, parameter_vector)  # MultinomialCGF$K1(tvec, c(0.2,0.3,0.5))
#> [1] 2 3 5
```
