#' Create a Poisson CGF Object
#'
#' @description
#' Creates a cumulant generating function (CGF) object for a Poisson-distributed random variable
#' with rate parameter `lambda`.
#'
#' The Poisson random variable \(X\) with parameter \(\lambda\) has:
#' \[
#' K(t) = \lambda(e^t - 1),
#' \]
#'
#'
#' @param lambda A positive numeric scalar. The rate parameter of the Poisson distribution.
#'
#' @return A `CGF` object (as created by `createCGF()`) specialized for the Poisson distribution.
#' It provides:
#' - `K(tvec, params)` returns \(\lambda(e^tvec - 1)\).
#' - `K1(tvec, params)` returns \(\lambda e^tvec\).
#' - `K2(tvec, params)` returns a matrix
#' - `K3operator(tvec, v1, v2, v3, params)` 
#' - `K4operator(tvec, v1, v2, v3, v4, params)` 
#'
#' All other methods of the `CGF` class (like `tilting_exponent`, `neg_ll`, and `func_T`) can be used on this object.
#'
#'
#' @noRd
createPoissonCGF <- function() {
  
  # Adaptor: extracts lambda from parameter_vector
  adaptor <- function(parameter_vector) {
    lambda <- parameter_vector[1]
    if (!is.numeric(lambda) || lambda <= 0) {
      stop("parameter_vector must be a numeric vector of length 1 with a positive value for lambda.")
    }
    lambda
  }
  
  K_vectorized <- function(tvec, lambda) { lambda * (exp(tvec) - 1) }
  K_func <- function(tvec, lambda) {   sum(K_vectorized(tvec, lambda)) }
  K1_func <- function(tvec, lambda) {  lambda * exp(tvec) }
  K2_vectorized <- function(tvec, lambda) { lambda * exp(tvec) }
  K2_func <- function(tvec, lambda) { diag(lambda * exp(tvec), nrow = length(tvec)) }
  
  K3_vectorized <- function(tvec, lambda) { lambda * exp(tvec) }
  K3operator_func <- function(tvec, v1, v2, v3, lambda) { sum(K3_vectorized(tvec, lambda) * v1 * v2 * v3) }
  
  K4_vectorized <- function(tvec, lambda) { lambda * exp(tvec) }
  K4operator_func <- function(tvec, v1, v2, v3, v4, lambda) {  sum(K4_vectorized(tvec, lambda) * v1 * v2 * v3 * v4)  }
  
  tilting_exponent_vectorized <- function(tvec, lambda) { K_vectorized(tvec, lambda) - tvec * K1(tvec, lambda) }
  tilting_exponent_func <- function(tvec, lambda) { sum(tilting_exponent_vectorized(tvec, lambda)) }
  
  neg_ll_func <- function(tvec, lambda) {
    sum(0.5 * log(2*pi*K2_vectorized(tvec, lambda)) - tilting_exponent_vectorized(tvec, lambda))
  }
  
  func_T_func <- function(tvec, lambda) {
    k2val <- K2_vectorized(tvec, lambda)
    k2sq_val <- k2val^2
    k3val <- K3_vectorized(tvec, lambda)
    k4val <- K4_vectorized(tvec, lambda)
    sum(k4val / (8 * k2sq_val) - 5 * k3val^2 / (24 * k2sq_val * k2val))
  }
  
  K4operatorAABB_func <- function(tvec, Q1, Q2, lambda) {
    sum(K4_vectorized(tvec, lambda) * diag(Q1) * diag(Q2))
  }
  K3K3operatorAABBCC_func <- function(tvec, Q1, Q2, Q3, lambda) {
    k3_vals <- K3_vectorized(tvec, lambda)
    sum( (diag(Q1) * k3_vals) %*% Q2 %*% (diag(Q3) * k3_vals) )
  }
  K3K3operatorABCABC_func <- function(tvec, Q1, Q2, Q3, lambda) {
    k3_vals <- K3_vectorized(tvec, lambda)
    mat_k3_vals <- diag(k3_vals, nrow = length(tvec))
    sum(mat_k3_vals %*% (Q1 %*% Q2 %*% Q3) %*% mat_k3_vals)
  }
  ineq_constraint_vectorized <- function(tvec, lambda) { numeric(0) }
  
  
  
  # Use createCGF to build the CGF object
  createCGF(
    K = K_func, 
    K1 = K1_func, 
    K2 = K2_func, 
    K3operator = K3operator_func, 
    K4operator = K4operator_func, 
    ineq_constraint = NULL,
    param_adaptor = adaptor,
    tilting_exponent = tilting_exponent_func,
    neg_ll = neg_ll_func,
    func_T = func_T_func,
    K2operator = NULL,
    K2operatorAK2AT = NULL,
    K4operatorAABB = K4operatorAABB_func,
    K3K3operatorAABBCC = K3K3operatorAABBCC_func,
    K3K3operatorABCABC = K3K3operatorABCABC_func,
    K4operatorAABB_factored = NULL,
    K3K3operatorAABBCC_factored = NULL,
    K3K3operatorABCABC_factored = NULL
  )
}
