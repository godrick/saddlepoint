# --------------------------------------------------------------------
# File: A-CGF_base.R
#
# PURPOSE & INSTRUCTIONS:
#  1) This file defines a base CGF (cumulant generating function) class using R6.
#  2) There are five compulsory methods any CGF must provide:
#       - K(tvec, parameter_vector)
#       - K1(tvec, parameter_vector)
#       - K2(tvec, parameter_vector)
#       - K3operator(tvec, parameter_vector, v1, v2, v3)
#       - K4operator(tvec, parameter_vector, v1, v2, v3, v4)
#  3) The class also supports optional methods (tilting_exponent, neg_ll, func_T, etc.). Users
#     can supply them or rely on defaults.
#  4) Some methods are private (e.g., neg_ll, func_T, and operator "factored" forms). We do not want 
#     them directly accessible via $ from an instance. They remain hidden in the private environment.
#  5) Public methods include a few optional operators (like K2operator) and the `.get_private_method`
#     for controlled access to certain private methods.
#  6) The `createCGF()` factory function at the end makes it easy to instantiate a CGF without subclassing.
#
#  NOTE ON ENVIRONMENTS:
#  - All names you intend to overwrite in the constructor must be declared
#    in 'public = list(...)' or 'private = list(...)' to avoid "cannot add bindings
#    to a locked environment" errors.
#  - By default, this class is not aware of the "param_adaptor". If you need parameter adaptation, see adaptCGF function.
#
#  EXTENSIBILITY:
#  - The 'initialize' method accepts '...' for additional named methods or overrides.
#  - End-users typically won't instantiate this class directly if they're using a specialized CGF,
#    but they can use createCGF(...) as a quick factory approach.
# --------------------------------------------------------------------




# --------------------------------------------------------------------
# Base CGF Class (R6)
# This is the base class for implementing CGF objects for various distributions.
# --------------------------------------------------------------------
#
# Required methods (subclass or createCGF construction):
#   1) K(tvec, parameter_vector)
#   2) K1(tvec, parameter_vector)
#   3) K2(tvec, parameter_vector)
#   4) K3operator(tvec, parameter_vector, v1, v2, v3)
#   5) K4operator(tvec, parameter_vector, v1, v2, v3, v4)
#
# Optional methods:
#   - tilting_exponent(tvec, parameter_vector)
#   - neg_ll(tvec, parameter_vector)
#   - func_T(tvec, parameter_vector)
#   - K2operator(tvec, parameter_vector, x, y)
#   - K2operatorAK2AT(tvec, parameter_vector, A)
#   - K4operatorAABB(tvec, parameter_vector, Q1, Q2)
#   - K3K3operatorAABBCC(tvec, parameter_vector, Q1, Q2, Q3)
#   - K3K3operatorABCABC(tvec, parameter_vector, Q1, Q2, Q3)
#   - K4operatorAABB_factored(tvec, parameter_vector, A1, d1, A2, d2)
#   - K3K3operatorAABBCC_factored(tvec, parameter_vector, A1, d1, A2, d2, A3, d3)
#   - K3K3operatorABCABC_factored(tvec, parameter_vector, A1, d1, A2, d2, A3, d3)
#   - ineq_constraint(tvec, parameter_vector)
#   - compute_analytic_tvec_hat(x, parameter_vector)
#
# Other objects:
#   - analytic_tvec_hat_func: a function that computes tvec from y and parameters
#   - op_name: a label for debugging/call history
#
# The 'initialize' method either uses the user-supplied or default implementation
# for each optional method.
# --------------------------------------------------------------------










# ------------------------------------------------------------------------
# Helper function: check_fun_sig
# This check that the user-supplied function has the right signature
# Currently not used.
# ------------------------------------------------------------------------
check_fun_sig <- function(fn, expected_args) {
  # If 'fn' is not a function, fail immediately
  if (!is.function(fn)) stop("Supplied object is not a function.")
  formals_list <- formals(fn)
  actual_args <- names(formals_list)
  # Compare argument names
  if (!identical(actual_args, expected_args)) {
    stop(
      "Invalid function signature. Expected arguments: ",
      paste(expected_args, collapse = ", "),
      "; got: ",
      paste(actual_args, collapse = ", ")
    )
  }
  # If we get here, 'fn' has the exact signature we want
}










#' @noRd
CGF <- R6::R6Class(
  classname = "CGF",
  
  # ----------------------------------------------------------
  # Private fields / methods
  # ----------------------------------------------------------
  private = list(
    
    # --- CORE user-supplied method pointers ---
    K_func = NULL,
    K1_func = NULL,
    K2_func = NULL,
    K3operator_func = NULL,
    K4operator_func = NULL,
    
    # --- OPTIONAL user-supplied method pointers ---
    ineq_constraint_func = NULL,
    analytic_tvec_hat_func = NULL,
    
    
    # "Hidden" or private-labeled methods:
    #  These are the default or user-supplied tilting_exponent, neg_ll, func_T, etc.
    tilting_exponent = NULL,
    neg_ll = NULL,
    func_T = NULL,
    K4operatorAABB_factored = NULL,
    K3K3operatorAABBCC_factored = NULL,
    K3K3operatorABCABC_factored = NULL,
    

    # Additional optional operator pointers
    K4operatorAABB_func = NULL,
    K3K3operatorAABBCC_func = NULL,
    K3K3operatorABCABC_func = NULL,
    K2operator_func = NULL,
    K2operatorAK2AT_func = NULL
  ),
  
  # ----------------------------------------------------------
  # Public members (accessible via $ on CGF object)
  # ----------------------------------------------------------
  public = list(
    
    # Keep a call_history for debugging/tracking
    call_history = NULL,
    
    # Pre-declare optional public methods here so we can safely overwrite them in `initialize`.
    K2operator = NULL,
    K2operatorAK2AT = NULL,
    K4operatorAABB  = NULL,
    K3K3operatorAABBCC = NULL,
    K3K3operatorABCABC = NULL,
    ineq_constraint = NULL,
    has_analytic_tvec_hat = NULL,
    analytic_tvec_hat = NULL,
    
    # -----------------------------------------------------------------------
    # CONSTRUCTOR
    # -----------------------------------------------------------------------
    initialize = function(K_func, K1_func, K2_func, K3operator_func, K4operator_func,
                          ineq_constraint_func = NULL,
                          # param_adaptor = function(x) x,
                          analytic_tvec_hat_func = NULL,
                          tilting_exponent_func = NULL,
                          neg_ll_func = NULL,
                          func_T_func = NULL,
                          K4operatorAABB_func = NULL,
                          K3K3operatorAABBCC_func = NULL,
                          K3K3operatorABCABC_func = NULL,
                          K4operatorAABB_factored_func = NULL,
                          K3K3operatorAABBCC_factored_func = NULL,
                          K3K3operatorABCABC_factored_func = NULL,
                          K2operator_func = NULL,
                          K2operatorAK2AT_func = NULL,
                          op_name = "UnnamedOperation",
                          ...
    ) {
      # --- Store core user-supplied methods ---
      private$K_func <- K_func
      private$K1_func <- K1_func
      private$K2_func <- K2_func
      private$K3operator_func <- K3operator_func
      private$K4operator_func <- K4operator_func
      
      
      # --- Store optional user-supplied methods ---
      private$ineq_constraint_func <- ineq_constraint_func
      private$analytic_tvec_hat_func <- analytic_tvec_hat_func
      private$K2operator_func <- K2operator_func
      private$K2operatorAK2AT_func <- K2operatorAK2AT_func
      private$K4operatorAABB_func <- K4operatorAABB_func
      private$K3K3operatorAABBCC_func <- K3K3operatorAABBCC_func
      private$K3K3operatorABCABC_func <- K3K3operatorABCABC_func
      
      
      
      # --- Assign or default for "factored" private methods ---
      if (!is.null(K4operatorAABB_factored_func)) {
        ##### We may use check_fun_sig here to ensure the user's function has the correct number of arguments
        ###   check_fun_sig(fn = K4operatorAABB_factored_func, expected_args = c("tvec", "parameter_vector", "A1", "d1", "A2", "d2"))
        ###   private$K4operatorAABB_factored <- K4operatorAABB_factored_func
        ##### But that check looks a bit strict.
        ##### For now we will simply wrap the supplied function which for now will ensure a universal signature.
        private$K4operatorAABB_factored <- function(tvec, parameter_vector, A1, d1, A2, d2) K4operatorAABB_factored_func(tvec, parameter_vector, A1, d1, A2, d2)
      } else {
        private$K4operatorAABB_factored <- function(tvec, parameter_vector, A1, d1, A2, d2) {
          r1 <- length(d1)
          r2 <- length(d2)
          res <- 0
          for (m1 in seq_len(r1)) {
            for (m2 in seq_len(r2)) {
              res <- res + d1[m1]*d2[m2]*self$K4operator(
                tvec, parameter_vector, A1[,m1], A1[,m1], A2[,m2], A2[,m2]
              )
            }
          }
          res
        }
      }
      
      if (!is.null(K3K3operatorAABBCC_factored_func)) {
        private$K3K3operatorAABBCC_factored <- function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) K3K3operatorAABBCC_factored_func(tvec, parameter_vector, A1, d1, A2, d2, A3, d3)
      } else {
        private$K3K3operatorAABBCC_factored <- function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
          r1 <- length(d1)
          r2 <- length(d2)
          r3 <- length(d3)
          res <- 0
          for (m2 in seq_len(r2)) {
            factor1 <- 0
            for (m1 in seq_len(r1)) {
              factor1 <- factor1 + d1[m1]*self$K3operator(tvec, parameter_vector, A1[,m1], A1[,m1], A2[,m2])
            }
            factor2 <- 0
            for (m3 in seq_len(r3)) {
              factor2 <- factor2 + d3[m3]*self$K3operator(tvec, parameter_vector, A2[,m2], A3[,m3], A3[,m3])
            }
            res <- res + d2[m2]*factor1*factor2
          }
          res
        }
      }
      
      if (!is.null(K3K3operatorABCABC_factored_func)) {
        private$K3K3operatorABCABC_factored <- function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) K3K3operatorABCABC_factored_func(tvec, parameter_vector, A1, d1, A2, d2, A3, d3)
      } else {
        private$K3K3operatorABCABC_factored <- function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
          r1 <- length(d1)
          r2 <- length(d2)
          r3 <- length(d3)
          message("The discrepancy option/compute.funcT has initiated a computation that may take a few moments...")
          res <- 0
          for (m1 in seq_len(r1)) {
            for (m2 in seq_len(r2)) {
              for (m3 in seq_len(r3)) {
                val <- self$K3operator(tvec, parameter_vector, A1[,m1], A2[,m2], A3[,m3])
                res <- res + d1[m1]*d2[m2]*d3[m3]*(val*val)
              }
            }
          }
          res
        }
      }
      
      # --- Assign or default for private tilting_exponent, neg_ll, and func_T ---
      if (!is.null(tilting_exponent_func)) {
        private$tilting_exponent <- function(tvec, parameter_vector) tilting_exponent_func(tvec, parameter_vector)
      } else {
        private$tilting_exponent <- function(tvec, parameter_vector) {
          self$K(tvec, parameter_vector) - sum(tvec * self$K1(tvec, parameter_vector))
        }
      }
      
      if (!is.null(neg_ll_func)) {
        private$neg_ll <- function(tvec, parameter_vector) neg_ll_func(tvec, parameter_vector)
      } else {
        private$neg_ll <- function(tvec, parameter_vector) {
          te <- private$tilting_exponent(tvec, parameter_vector)
          K2_val <- self$K2(tvec, parameter_vector)
          val_logdet <- determinant(K2_val, logarithm = TRUE)$modulus
          0.5 * val_logdet + 0.5 * length(tvec)*log(2*pi) - te
        }
      }

      if (!is.null(func_T_func)) {
        private$func_T <- function(tvec, parameter_vector) func_T_func(tvec, parameter_vector)
      } else {
        private$func_T <- function(tvec, parameter_vector) {
          K2_val <- self$K2(tvec, parameter_vector)
          K2_inv <- solve(K2_val)
          chol_K2_inv <- chol(K2_inv)
          diag_K2_inv <- diag(chol_K2_inv)
          d <- diag_K2_inv^2
          A <- t(chol_K2_inv) %*% diag(1/diag_K2_inv)
          
          K4_AABB   <- private$K4operatorAABB_factored(tvec, parameter_vector, A, d, A, d)
          K3K3_ABBC <- private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A, d, A, d, A, d)
          K3K3_ABC  <- private$K3K3operatorABCABC_factored(tvec, parameter_vector, A, d, A, d, A, d)
          K4_AABB/8 - K3K3_ABBC/8 - K3K3_ABC/12
        }
      }

      
      # --- Record operation name in the call history ---
      if (!is.character(op_name)) stop("'operation' must be of type character")
      self$call_history <- if (!is.null(self$call_history)) {
        c(self$call_history, op_name)
      } else {
        op_name
      }
      
      
      # ---------------------------------------------------------------------
      # Overwrite the public optional operators if the user supplied custom versions
      # Otherwise, assign defaults
      #---------------------------------------------------------------------
      
      # K2operator
      if (!is.null(K2operator_func)) {
        self$K2operator <- function(tvec, parameter_vector, x, y) {
          K2operator_func(tvec, parameter_vector, x, y)
        }
      } else {
        self$K2operator <- function(tvec, parameter_vector, x, y) {
          K2_val <- self$K2(tvec, parameter_vector)
          as.vector(t(x) %*% (K2_val %*% y))
        }
      }
      
      # K2operatorAK2AT
      if (!is.null(K2operatorAK2AT_func)) {
        self$K2operatorAK2AT <- function(tvec, parameter_vector, A) {
          K2operatorAK2AT_func(tvec, parameter_vector, A)
        }
      } else {
        self$K2operatorAK2AT <- function(tvec, parameter_vector, A) {
          K2_val <- self$K2(tvec, parameter_vector)
          A %*% K2_val %*% t(A)
        }
      }
      
      # K4operatorAABB
      if (!is.null(K4operatorAABB_func)) {
        self$K4operatorAABB <- function(tvec, parameter_vector, Q1, Q2) {
          K4operatorAABB_func(tvec, parameter_vector, Q1, Q2)
        }
      } else {
        self$K4operatorAABB <- function(tvec, parameter_vector, Q1, Q2) {
          chol_Q1 <- chol(Q1)
          diag_Q1 <- diag(chol_Q1)
          d1 <- diag_Q1^2
          A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
          private$K4operatorAABB_factored(tvec, parameter_vector, A1, d1, A1, d1) 
        }
      }
      
      # K3K3operatorAABBCC
      if (!is.null(K3K3operatorAABBCC_func)) {
        self$K3K3operatorAABBCC <- function(tvec, parameter_vector, Q1, Q2, Q3) {
          K3K3operatorAABBCC_func(tvec, parameter_vector, Q1, Q2, Q3)
        }
      } else {
        self$K3K3operatorAABBCC <- function(tvec, parameter_vector, Q1, Q2, Q3) {
          chol_Q1 <- chol(Q1)
          diag_Q1 <- diag(chol_Q1)
          d1 <- diag_Q1^2
          A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
          private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1)
        }
      }
      
      # K3K3operatorABCABC
      if (!is.null(K3K3operatorABCABC_func)) {
        self$K3K3operatorABCABC <- function(tvec, parameter_vector, Q1, Q2, Q3) {
          K3K3operatorABCABC_func(tvec, parameter_vector, Q1, Q2, Q3)
        }
      } else {
        self$K3K3operatorABCABC <- function(tvec, parameter_vector, Q1, Q2, Q3) {
          chol_Q1 <- chol(Q1)
          diag_Q1 <- diag(chol_Q1)
          d1 <- diag_Q1^2
          A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
          private$K3K3operatorABCABC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1)
        }
      }
      
      # ineq_constraint
      if (!is.null(ineq_constraint_func)) {
        self$ineq_constraint <- function(tvec, parameter_vector) {
          ineq_constraint_func(tvec, parameter_vector)
        }
      } else {
        self$ineq_constraint <- function(tvec, parameter_vector) {
          numeric(0)
        }
      }
      
      # has_analytic_tvec_hat / analytic_tvec_hat
      self$has_analytic_tvec_hat <- function() {
        !is.null(private$analytic_tvec_hat_func)
      }
      if(self$has_analytic_tvec_hat()) {
        self$analytic_tvec_hat <- function(x, parameter_vector) {
          stopifnot(is.numeric(x), !any(x <= 0))
          private$analytic_tvec_hat_func(x, parameter_vector)
        }
      } else {
        self$analytic_tvec_hat <- NULL
      }
        
    },
    
    
    
    
    # -----------------------------------------------------------------------
    # REQUIRED PUBLIC METHODS 
    # -----------------------------------------------------------------------
    K = function(tvec, parameter_vector) {
      private$K_func(tvec, parameter_vector)
    },
    
    K1 = function(tvec, parameter_vector) {
      private$K1_func(tvec, parameter_vector)
    },
    
    K2 = function(tvec, parameter_vector) {
      private$K2_func(tvec, parameter_vector)
    },
    
    K3operator = function(tvec, parameter_vector, v1, v2, v3) {
      private$K3operator_func(tvec, parameter_vector, v1, v2, v3)
    },
    
    K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) {
      private$K4operator_func(tvec, parameter_vector, v1, v2, v3, v4)
    },
    
    
    
    
    
    # -----------------------------------------------------------------------
    # Print method
    # -----------------------------------------------------------------------
    print = function(...) {
      cat("<CGF Object>\n")
      if (!is.null(self$call_history)) {
        mapping_str <- paste(self$call_history, collapse = " -> ")
        cat("Used:", mapping_str, "\n")
      }
      cat("Class hierarchy:", paste(class(self), collapse = " -> "), "\n")
      invisible(self)
    },
    
    # -----------------------------------------------------------------------
    # CONTROLLED ACCESS to Private Methods
    # -----------------------------------------------------------------------
    .get_private_method = function(method_name) {
      if (!is.character(method_name) || length(method_name) != 1) {
        stop("'method_name' must be a single character string.")
      }
      if (!method_name %in% names(private)) {
        stop(paste0("'", method_name, "' is not a private method in the CGF class."))
      }
      # a whitelist of allowed private methods
      allowed_methods <- c(
        "neg_ll", 
        "tilting_exponent", 
        "K4operatorAABB_factored", 
        "K3K3operatorAABBCC_factored", 
        "K3K3operatorABCABC_factored",
        "func_T"
      )
      if (!(method_name %in% allowed_methods)) {
        stop(paste0("Access to private method '", method_name, "' is not permitted."))
      }
      
      method_ <- private[[method_name]]
      if (!is.function(method_)) {
        stop(paste0("Private member '", method_name, "' is not a function."))
      }
      method_
    },
    
    # -----------------------------------------------------------------------
    # EXAMPLE Additional Method: compute.spa.negll
    # (calls a global compute.spa.negll function with the cgf object)
    # -----------------------------------------------------------------------
    compute.spa.negll = function(parameter_vector,
                                 observed.data,
                                 tvec.hat = NULL,
                                 gradient = FALSE,
                                 hessian = FALSE,
                                 ...) {
      compute.spa.negll(
        cgf              = self,
        parameter_vector = parameter_vector,
        observed.data    = observed.data,
        tvec.hat         = tvec.hat,
        gradient         = gradient,
        hessian          = hessian,
        ...
      )
    }
  )
)



# ------------------------------------------------------------------------
# FACTORY FUNCTION: createCGF
# ------------------------------------------------------------------------
#' Create a CGF object from user-defined functions
#'
#' 
#' @description
#' This creates an object of type CGF using user-supplied functions. You supply 
#' the five essential methods (`K`, `K1`, `K2`, `K3operator`, `K4operator`) plus 
#' any optional overrides (e.g., `tilting_exponent` or `neg_ll`), and it returns
#' a `CGF` instance. 
#'
#'
#' @param K A function `K(tvec, parameter_vector) -> numeric scalar`.
#' @param K1 A function `K1(tvec, parameter_vector) -> numeric vector`.
#' @param K2 A function `K2(tvec, parameter_vector) -> numeric matrix`.
#' @param K3operator A function implementing the third-order operator.
#' @param K4operator A function implementing the fourth-order operator.
#' 
#' @param ineq_constraint Optional function for inequality constraints.
#' 
#' @param analytic_tvec_hat_func Optional function for an analytic solution
#'   of the saddlepoint equation. If provided, call it via `cgf$analytic_tvec_hat(x, param)`.
#' @param op_name A descriptive label for this CGF object. Used for the print method.
#'               Default is "UnnamedOperation".
#' 
#' @param tilting_exponent (optional) Overriding function for the tilting exponent.
#' @param neg_ll (optional) Overriding function for the negative log-likelihood.
#' @param func_T (optional) Overriding function for the first-order correction term.
#' @param K2operator,K2operatorAK2AT,K4operatorAABB,K3K3operatorAABBCC,K3K3operatorABCABC (optional) Overriding operator methods.
#' @param K4operatorAABB_factored,K3K3operatorAABBCC_factored,K3K3operatorABCABC_factored (optional) Overriding factored-operator methods.
#' @param ... Additional named methods or overrides. This demonstrates the extensibility
#'   mechanism, letting you add new CGF-related methods without altering the base class.              
#'
#' @return An object of class `CGF`.
#' @export
createCGF <- function(K, K1, K2, K3operator, K4operator, 
                      ineq_constraint = NULL,
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
  # Collect optional methods into a list
  user_optional_methods <- list(
    tilting_exponent_func          = tilting_exponent,
    neg_ll_func                    = neg_ll,
    func_T_func                    = func_T,
    K4operatorAABB_func           = K4operatorAABB,
    K3K3operatorAABBCC_func       = K3K3operatorAABBCC,
    K3K3operatorABCABC_func       = K3K3operatorABCABC,
    K4operatorAABB_factored_func  = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored_func = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored_func = K3K3operatorABCABC_factored,
    K2operator_func               = K2operator,
    K2operatorAK2AT_func          = K2operatorAK2AT
  )
  
  # Any additional named overrides passed via ...
  additional_methods <- list(...)
  
  # Merge user-supplied with additional; latter has precedence
  all_optional_methods <- modifyList(user_optional_methods, additional_methods)
  
  # Construct and return the CGF object
  do.call(CGF$new, c(
    list(
      K_func               = K,
      K1_func              = K1,
      K2_func              = K2,
      K3operator_func      = K3operator,
      K4operator_func      = K4operator,
      ineq_constraint_func = ineq_constraint,
      analytic_tvec_hat_func = analytic_tvec_hat_func,
      op_name = op_name
    ),
    all_optional_methods
  ))
}
