# Create a `CGF` object from vectorized functions

This function allows you to create a `CGF` object using vectorized
functions along with any optional operators or methods.

## Usage

``` r
createCGF_fromVectorisedFunctions(
  K_vectorized_func,
  K1_vectorized_func,
  K2_vectorized_func,
  K3_vectorized_func,
  K4_vectorized_func,
  ineq_constraint = NULL,
  analytic_tvec_hat = NULL,
  op_name = "UnnamedOperation",
  tilting_exponent = NULL,
  neg_ll = NULL,
  func_T = NULL,
  rsim = NULL,
  K4operatorAABB = NULL,
  K3K3operatorAABBCC = NULL,
  K3K3operatorABCABC = NULL,
  K4operatorAABB_factored = NULL,
  K3K3operatorAABBCC_factored = NULL,
  K3K3operatorABCABC_factored = NULL,
  K2operator = NULL,
  K2operatorAK2AT = NULL,
  ...
)
```

## Arguments

- K_vectorized_func:

  A function of the form `function(tvec, param) -> numeric vector` that
  returns the CGF values.

- K1_vectorized_func:

  A function of the form `function(tvec, param) -> numeric vector` that
  returns the first derivative values.

- K2_vectorized_func:

  A function of the form `function(tvec, param) -> numeric vector` that
  returns the second derivative values.

- K3_vectorized_func, K4_vectorized_func:

  Similar vectorized functions for the third and fourth derivatives.

- ineq_constraint:

  Optional inequality constraint function.

- analytic_tvec_hat:

  Optional `tvec_hat` function override.

- op_name:

  Optional character string indicating the name of the operation or
  transformation being performed.

- tilting_exponent:

  Optional tilting exponent function override.

- neg_ll:

  Optional neg_ll function override.

- func_T:

  Optional func_T function override.

- rsim:

  Optional simulation method. A function of the form
  `function(n, vector_length, parameter_vector, tvec = NULL, ...)`
  returning a numeric vector or matrix of length `n * vector_length` or
  a `vector_length x n` matrix.

- K4operatorAABB, K3K3operatorAABBCC, K3K3operatorABCABC:

  Optional operator overrides.

- K4operatorAABB_factored, K3K3operatorAABBCC_factored,
  K3K3operatorABCABC_factored:

  Optional factored operator overrides.

- K2operator, K2operatorAK2AT:

  Optional operator overrides.

- ...:

  Any additional named optional methods.

## Value

A `CGF` object.
