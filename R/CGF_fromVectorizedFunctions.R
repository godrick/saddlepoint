# R/CGF_fromVectorizedFunctions.R

#' A Vectorized CGF Class
#'
#' @description
#' `VectorOfIIDCGF` inherits from `CGF` and simplifies creation of a CGF object
#' when you have functions that have vectorized forms.
#' 
#' By providing vectorized versions of K, K1, K2, K3, K4 (and possibly their operators),
#' `VectorOfIIDCGF` can construct a standard `CGF` object.
#'
#' For instance:
#' - `K_vectorized(tvec, params)` returns a vector, one entry per t in tvec.
#'   `K()` will be defined as `sum(K_vectorized(tvec, params))`.
#'
#' - `K1_vectorized(tvec, params)` returns a vector the same length as tvec (the gradient).
#'   `K1()` can directly return this result.
#'
#' - `K2_vectorized(tvec, params)` returns a vector (one per t) that represents the diagonal
#'   elements of K2 if it's diagonal. If not, you can provide more complex logic.
#'
#' Similar logic applies for K3 and K4 operator forms.
#'
#' @noRd
VectorOfIIDCGF <- R6::R6Class("VectorOfIIDCGF",
                             inherit = CGF,
                             
                             public = list(
                               K_vectorized_func = NULL,
                               K1_vectorized_func = NULL,
                               K2_vectorized_func = NULL,
                               K3_vectorized_func = NULL,
                               K4_vectorized_func = NULL,
                               
                               initialize = function(
                                      K_vectorized_func, K1_vectorized_func, K2_vectorized_func,
                                      K3_vectorized_func, K4_vectorized_func,
                                      ineq_constraint_func = NULL,
                                      param_adaptor = function(x) x,
                                      ...
                               ) {
                                 # Store vectorized functions
                                 self$K_vectorized_func <- K_vectorized_func
                                 self$K1_vectorized_func <- K1_vectorized_func
                                 self$K2_vectorized_func <- K2_vectorized_func
                                 self$K3_vectorized_func <- K3_vectorized_func
                                 self$K4_vectorized_func <- K4_vectorized_func
                                 
                                 # Define base K, K1, K2, K3operator, K4operator from vectorized forms if not provided
                                 K_func <- function(tvec, p) {
                                   # sum of K_vectorized over tvec
                                   sum(self$K_vectorized_func(tvec, p))
                                 }
                                 K1_func <- function(tvec, p) {
                                   # K1 is already a vector over tvec
                                   self$K1_vectorized_func(tvec, p)
                                 }
                                 K2_func <- function(tvec, p) {
                                   # If K2 is diagonal for IID distributions, we can form a diag matrix:
                                   diag(self$K2_vectorized_func(tvec, p), nrow = length(tvec))
                                 }
                                 
                                 
                                 # sum(K3_vectorized(tvec)*v1*v2*v3)
                                 K3operator_func <- function(tvec, p, v1, v2, v3) {
                                   vals <- self$K3_vectorized_func(tvec, p)
                                   sum(vals * v1 * v2 * v3)
                                 }
                                 
                                 
                                 # sum(K4_vectorized(tvec)*v1*v2*v3*v4)
                                 K4operator_func <- function(tvec, p, v1, v2, v3, v4) {
                                   vals <- self$K4_vectorized_func(tvec, p)
                                   sum(vals * v1 * v2 * v3 * v4)
                                 }
                                 
                                 tilting_exponent_func <- function(tvec, p) {
                                   K_vals <- self$K_vectorized_func(tvec, p)
                                   K1_vals <- self$K1_vectorized_func(tvec, p)
                                   sum(K_vals - tvec * K1_vals)
                                 }
                                 
                                 
                                 neg_ll_func <- function(tvec, p) {
                                   K2_vals <- self$K2_vectorized_func(tvec, p)
                                   K_vals <- self$K_vectorized_func(tvec, p)
                                   K1_vals <- self$K1_vectorized_func(tvec, p)
                                   tilting_vals <- K_vals - tvec * K1_vals
                                   sum(0.5 * log(2*pi*K2_vals) - tilting_vals)
                                 }
                                 
                                 func_T_func <- function(tvec, p) {
                                   k2val <- self$K2_vectorized_func(tvec, p)
                                   k2sq_val <- k2val^2
                                   k3val <- self$K3_vectorized_func(tvec, p)
                                   k4val <- self$K4_vectorized_func(tvec, p)
                                   sum(k4val / (8 * k2sq_val) - 5 * k3val^2 / (24 * k2sq_val * k2val))
                                 }
                                 
                                 K4operatorAABB_func <- function(tvec, Q1, Q2, p) {
                                   sum(self$K4_vectorized_func(tvec, p) * diag(Q1) * diag(Q2))
                                 }
                                 
                                 K3K3operatorAABBCC_func <- function(tvec, Q1, Q2, Q3, p) {
                                   k3_vals <- self$K3_vectorized_func(tvec, p)
                                   sum( (diag(Q1) * k3_vals) %*% Q2 %*% (diag(Q3) * k3_vals) )
                                 }
                                 
                                 K3K3operatorABCABC_func <- function(tvec, Q1, Q2, Q3, p) {
                                   k3_vals <- self$K3_vectorized_func(tvec, p)
                                   mat_k3_vals <- diag(k3_vals, nrow = length(tvec))
                                   sum(mat_k3_vals %*% (Q1 %*% Q2 %*% Q3) %*% mat_k3_vals)
                                 }
                                 
                                 # Call super$initialize() with these derived functions
                                 super$initialize(
                                   K_func = K_func,
                                   K1_func = K1_func,
                                   K2_func = K2_func,
                                   K3operator_func = K3operator_func,
                                   K4operator_func = K4operator_func,
                                   ineq_constraint_func = ineq_constraint_func,
                                   param_adaptor = param_adaptor,
                                   tilting_exponent_func = tilting_exponent_func,
                                   neg_ll_func = neg_ll_func,
                                   func_T_func = func_T_func,
                                   K4operatorAABB_func = K4operatorAABB_func,
                                   K3K3operatorAABBCC_func = K3K3operatorAABBCC_func,
                                   K3K3operatorABCABC_func = K3K3operatorABCABC_func,
                                   ...
                                 )
                               }
                             )
)
