# R/CGF_fromVectorizedFunctions.R

#' A Vectorized CGF Class
#'
#' @description
#' `VectorizedCGF` inherits from `CGF` and simplifies creation of a CGF object
#' when you have functions that have vectorized forms.
#' 
#' By providing vectorized versions of K, K1, K2, K3, K4 (and possibly their operators),
#' `VectorizedCGF` can construct a standard `CGF` object.
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
#' @export
VectorizedCGF <- R6::R6Class("VectorizedCGF",
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
    param_adaptor = function(x) x,
    # Optional overrides:
    K3operator_func = NULL, K4operator_func = NULL,
    tilting_exponent_func = NULL, neg_ll_func = NULL, func_T_func = NULL,
    K2operator_func = NULL, K2operatorAK2AT_func = NULL,
    K4operatorAABB_func = NULL, K3K3operatorAABBCC_func = NULL, K3K3operatorABCABC_func = NULL,
    K4operatorAABB_factored_func = NULL, K3K3operatorAABBCC_factored_func = NULL, K3K3operatorABCABC_factored_func = NULL,
    ineq_constraint_func = NULL
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
                                 
                                 # Default K3operator if not provided:
                                 # K3(tvec,v1,v2,v3) = sum(K3_vectorized(tvec)*v1*v2*v3)
                                 if (is.null(K3operator_func)) {
                                   K3operator_func <- function(tvec, p, v1, v2, v3) {
                                     vals <- self$K3_vectorized_func(tvec, p)
                                     sum(vals * v1 * v2 * v3)
                                   }
                                 }
                                 
                                 # Default K4operator if not provided:
                                 # K4(tvec,v1,v2,v3,v4) = sum(K4_vectorized(tvec)*v1*v2*v3*v4)
                                 if (is.null(K4operator_func)) {
                                   K4operator_func <- function(tvec, p, v1, v2, v3, v4) {
                                     vals <- self$K4_vectorized_func(tvec, p)
                                     sum(vals * v1 * v2 * v3 * v4)
                                   }
                                 }
                                 
                                 # Call super$initialize() with these derived functions
                                 super$initialize(
                                   K_func = K_func,
                                   K1_func = K1_func,
                                   K2_func = K2_func,
                                   K3operator_func = K3operator_func,
                                   K4operator_func = K4operator_func,
                                   ineq_constraint_func = ineq_constraint_func,
                                   tilting_exponent_func = tilting_exponent_func,
                                   neg_ll_func = neg_ll_func,
                                   func_T_func = func_T_func,
                                   K2operator_func = K2operator_func,
                                   K2operatorAK2AT_func = K2operatorAK2AT_func,
                                   K4operatorAABB_func = K4operatorAABB_func,
                                   K3K3operatorAABBCC_func = K3K3operatorAABBCC_func,
                                   K3K3operatorABCABC_func = K3K3operatorABCABC_func,
                                   K4operatorAABB_factored_func = K4operatorAABB_factored_func,
                                   K3K3operatorAABBCC_factored_func = K3K3operatorAABBCC_factored_func,
                                   K3K3operatorABCABC_factored_func = K3K3operatorABCABC_factored_func,
                                   param_adaptor = param_adaptor
                                 )
                               }
                             )
)
