# R/A-CGF_fromVectorizedFunctions.R
# Objects: VectorizedFunctionsCGF, createCGF_fromVectorisedFunctions

# The Setup of the followiing class.
# Compulsory Methods: 
# The child class VectorizedFunctionsCGF requires the five vectorized functions 
# (K_vectorized_func, K1_vectorized_func, etc.) as explicit arguments in its initialize method. 
#
#
# Child Defaults for Optional Methods:
# For each optional method (tilting_exponent_func, neg_ll_func, etc.), 
# the child class checks if the user has provided an override (i.e., if the argument is not NULL).
# If the user does not provide an override (NULL), the child class assigns its improved default implementation.
# If the user does provide an override (including NULL), it is passed to the parent class, 
# allowing users to either use their own function or revert to the parent’s default by setting it to NULL.
#
#
# Passing to Parent Class
# The child class constructs the required base methods (K_func, K1_func, etc.) from the vectorized functions.
# It then calls super$initialize() with all the required and optional methods.
# The ... mechanism allows for any additional optional methods to be seamlessly passed to the parent class, 
# supporting future extensibility.
# This ensures that new methods can be added in the future without modifying the child class.








#' A Vectorized CGF Class
#'
#' @description
#' `VectorizedFunctionsCGF` inherits from `CGF` and makes it easier to create a CGF
#' when you have vectorized functions.
#'
#' This class constructs the necessary `K()`, `K1()`, `K2()`, `K3operator()`, 
#' and `K4operator()` from your vectorized forms. It also sets default 
#' optional methods (tilting_exponent, neg_ll, func_T, etc.) based on 
#' the vectorized logic.
#'
#'
#' @noRd
VectorizedFunctionsCGF <- R6::R6Class(
  "VectorizedFunctionsCGF",
  inherit = CGF,  # base CGF
  
  private = list(
    K_vectorized_func  = NULL,
    K1_vectorized_func = NULL,
    K2_vectorized_func = NULL,
    K3_vectorized_func = NULL,
    K4_vectorized_func = NULL,
    
    # Internal helper method(s) to set child defaults for optional methods
    set_child_defaults = function(
              tilting_exponent_func,
              neg_ll_func,
              func_T_func,
              K4operatorAABB_func,
              K3K3operatorAABBCC_func,
              K3K3operatorABCABC_func
    ) {
      # If user didn't supply a function, create a default implementation
      if (is.null(tilting_exponent_func)) {
        tilting_exponent_func <- function(tvec, p) {
          K_vals <- private$K_vectorized_func(tvec, p)
          K1_vals <- private$K1_vectorized_func(tvec, p)
          sum(K_vals - tvec * K1_vals)
        }
      }
      if (is.null(neg_ll_func)) {
        neg_ll_func <- function(tvec, p) {
          K2_vals <- private$K2_vectorized_func(tvec, p)
          K_vals <- private$K_vectorized_func(tvec, p)
          K1_vals <- private$K1_vectorized_func(tvec, p)
          tilting_vals <- K_vals - tvec * K1_vals
          sum(0.5 * log(2*pi*K2_vals) - tilting_vals)
        }
      }
      if(is.null(func_T_func)){
        func_T_func <- function(tvec, p) {
          k2val <- private$K2_vectorized_func(tvec, p)
          k2sq_val <- k2val^2
          k3val <- private$K3_vectorized_func(tvec, p)
          k4val <- private$K4_vectorized_func(tvec, p)
          sum(k4val / (8 * k2sq_val) - 5 * k3val^2 / (24 * k2sq_val * k2val))
        }
      }
      if(is.null(K4operatorAABB_func)){
        K4operatorAABB_func <- function(tvec, p, Q1, Q2) {
          sum(private$K4_vectorized_func(tvec, p) * diag(Q1) * diag(Q2))
        }
      }
      if(is.null(K3K3operatorAABBCC_func)){
        K3K3operatorAABBCC_func <- function(tvec, p, Q1, Q2, Q3) {
          k3_vals <- private$K3_vectorized_func(tvec, p)
          sum( (diag(Q1) * k3_vals) %*% Q2 %*% (diag(Q3) * k3_vals) )
        }
      }
      if(is.null(K3K3operatorABCABC_func)){
        K3K3operatorABCABC_func <- function(tvec, p, Q1, Q2, Q3) {
          k3_vals <- private$K3_vectorized_func(tvec, p)
          mat_k3_vals <- diag(k3_vals, nrow = length(tvec))
          sum(mat_k3_vals %*% (Q1 %*% Q2 %*% Q3) %*% mat_k3_vals)
        }
      }
      
      list(
        tilting_exponent_func = tilting_exponent_func,
        neg_ll_func           = neg_ll_func,
        func_T_func           = func_T_func,
        K4operatorAABB_func   = K4operatorAABB_func,
        K3K3operatorAABBCC_func = K3K3operatorAABBCC_func,
        K3K3operatorABCABC_func = K3K3operatorABCABC_func
      )
    }
  ),
  
  public = list(
    initialize = function(
                K_vectorized_func,
                K1_vectorized_func,
                K2_vectorized_func,
                K3_vectorized_func,
                K4_vectorized_func,
                
                # Optional overrides
                ineq_constraint_func = NULL,
                analytic_tvec_hat_func        = NULL,
                op_name = "UnnamedOperation",
                
                tilting_exponent_func = NULL,
                neg_ll_func           = NULL,
                func_T_func           = NULL,
                K4operatorAABB_func     = NULL,
                K3K3operatorAABBCC_func = NULL,
                K3K3operatorABCABC_func = NULL,
                K4operatorAABB_factored_func = NULL,
                K3K3operatorAABBCC_factored_func = NULL,
                K3K3operatorABCABC_factored_func = NULL,
                K2operator_func = NULL,
                K2operatorAK2AT_func = NULL,
                ...
    ) {
      # Store vectorized functions in private fields
      private$K_vectorized_func  <- K_vectorized_func
      private$K1_vectorized_func <- K1_vectorized_func
      private$K2_vectorized_func <- K2_vectorized_func
      private$K3_vectorized_func <- K3_vectorized_func
      private$K4_vectorized_func <- K4_vectorized_func
      
      # Call child default logic
      child_defaults <- private$set_child_defaults(
        tilting_exponent_func  = tilting_exponent_func,
        neg_ll_func            = neg_ll_func,
        func_T_func            = func_T_func,
        K4operatorAABB_func    = K4operatorAABB_func,
        K3K3operatorAABBCC_func= K3K3operatorAABBCC_func,
        K3K3operatorABCABC_func= K3K3operatorABCABC_func
      )
      
      
      super$initialize(
        K_func                = function(tvec, p) { sum(private$K_vectorized_func(tvec, p))  },
        K1_func               = function(tvec, p) { private$K1_vectorized_func(tvec, p)  },
        K2_func               = function(tvec, p) { diag(private$K2_vectorized_func(tvec, p), nrow = length(tvec)) },
        K3operator_func       = function(tvec, p, v1, v2, v3) { sum(private$K3_vectorized_func(tvec, p) * v1 * v2 * v3) },
        K4operator_func       = function(tvec, p, v1, v2, v3, v4) { sum(private$K4_vectorized_func(tvec, p) * v1 * v2 * v3 * v4) },
        
        ineq_constraint_func  = ineq_constraint_func,
        analytic_tvec_hat_func= analytic_tvec_hat_func,
        op_name        = op_name,
        
        # optional child-defaulted methods
        tilting_exponent_func = child_defaults$tilting_exponent_func,
        neg_ll_func           = child_defaults$neg_ll_func,
        func_T_func           = child_defaults$func_T_func,
        
        K4operatorAABB_func          = child_defaults$K4operatorAABB_func,
        K3K3operatorAABBCC_func      = child_defaults$K3K3operatorAABBCC_func,
        K3K3operatorABCABC_func      = child_defaults$K3K3operatorABCABC_func,
        
        K4operatorAABB_factored_func          = K4operatorAABB_factored_func,
        K3K3operatorAABBCC_factored_func      = K3K3operatorAABBCC_factored_func,
        K3K3operatorABCABC_factored_func      = K3K3operatorABCABC_factored_func,
        
        K2operator_func           = K2operator_func,
        K2operatorAK2AT_func      = K2operatorAK2AT_func,
        ...  # pass any further methods to the parent
      )
    }
    
    
    
  )
)




 
 







#' Create a `CGF` object from vectorized functions
#'
#' @description
#' This function allows you to create a `CGF` object using vectorized functions
#' along with any optional operators or methods.
#'
#' @param K_vectorized_func Function: returns a vector of K-values for tvec.
#' @param K1_vectorized_func Function: returns a vector (gradient) for K1.
#' @param K2_vectorized_func Function: returns a vector of second derivatives (often diagonal for IID).
#' @param K3_vectorized_func,K4_vectorized_func Similar vectorized functions for K3 and K4.
#' @param ineq_constraint_func Optional inequality constraint function.
#' @param analytic_tvec_hat_func Optional tvec_hat function override.
#' @param op_name Optional character string indicating the name of the operation or transformation being performed.
#' @param tilting_exponent Optional tilting exponent function override.
#' @param neg_ll Optional neg_ll function override.
#' @param func_T Optional func_T function override.
#' @param K4operatorAABB,K3K3operatorAABBCC,K3K3operatorABCABC Optional operator overrides.
#' @param K4operatorAABB_factored,K3K3operatorAABBCC_factored,K3K3operatorABCABC_factored Optional factored operator overrides.
#' @param K2operator,K2operatorAK2AT Optional operator overrides.
#' @param ... Any additional optional methods.
#'
#' @return A `CGF` object.
#' @export
createCGF_fromVectorisedFunctions <- function(
    K_vectorized_func,
    K1_vectorized_func,
    K2_vectorized_func,
    K3_vectorized_func,
    K4_vectorized_func,
    ineq_constraint_func = NULL,
    analytic_tvec_hat_func = NULL,
    op_name = "UnnamedOperation",
    tilting_exponent = NULL,
    neg_ll = NULL,
    func_T = NULL,
    K4operatorAABB = NULL,
    K3K3operatorAABBCC = NULL,
    K3K3operatorABCABC = NULL,
    K4operatorAABB_factored = NULL,
    K3K3operatorAABBCC_factored = NULL,
    K3K3operatorABCABC_factored = NULL,
    K2operator = NULL,
    K2operatorAK2AT = NULL,
    ...
) {
  
  # user-supplied optional methods 
  user_optional_methods <- list(
    tilting_exponent_func       = tilting_exponent,
    neg_ll_func                 = neg_ll,
    func_T_func                 = func_T,
    K4operatorAABB_func         = K4operatorAABB,
    K3K3operatorAABBCC_func     = K3K3operatorAABBCC,
    K3K3operatorABCABC_func     = K3K3operatorABCABC,
    K4operatorAABB_factored_func    = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored_func= K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored_func= K3K3operatorABCABC_factored,
    K2operator_func            = K2operator,
    K2operatorAK2AT_func       = K2operatorAK2AT
  )
  
  
  # any additional methods passed via ...
  additional_methods <- list(...)
  
  # Merge user-supplied optional methods with additional methods
  # Additional methods take precedence in case of name conflicts
  all_optional_methods <- modifyList(user_optional_methods, additional_methods)
  
  
  do.call(VectorizedFunctionsCGF$new, c(
    list(
      K_vectorized_func  = K_vectorized_func,
      K1_vectorized_func = K1_vectorized_func,
      K2_vectorized_func = K2_vectorized_func,
      K3_vectorized_func = K3_vectorized_func,
      K4_vectorized_func = K4_vectorized_func,
      ineq_constraint_func          = ineq_constraint_func,
      analytic_tvec_hat_func        = analytic_tvec_hat_func,
      op_name  = op_name
    ),
    all_optional_methods
  ))
  
}

