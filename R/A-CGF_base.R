# --------------------------------------------------------------------
# File: A-CGF_base.R
#
# Purpose: Provides a base CGF (cumulant generating function) class
#          supporting essential methods (K, K1, K2, K3operator, K4operator).
#          and optional ones (tilting_exponent, neg_ll, func_T, etc.).
#          Optional methods (such as tilting_exponent, neg_ll, and various operator forms)
#          can be supplied or will default to a basic implementation.
#          Also includes a factory function 'createCGF()' for convenience.
#
# Note: -This class is NOT exported. 
#        End-users will use specialized CGF objects or the createCGF() function.
#       -This class is not aware of parameter adaptation (param_adaptor is neutral/identity). If you want it, 
#        see the AdaptCGF class.
# 
#  
# --------------------------------------------------------------------



# --------------------------------------------------------------------
# Base CGF Class (R6)
# This is the base class for implementing CGF objects for various distributions.
# --------------------------------------------------------------------
#
# Required methods in any subclass or in the 'createCGF' construction:
#   1) K(tvec, param)
#   2) K1(tvec, param)
#   3) K2(tvec, param)
#   4) K3operator(tvec, param, v1, v2, v3)
#   5) K4operator(tvec, param, v1, v2, v3, v4)
#
# Optional methods may be overridden:
#   - tilting_exponent(tvec, param)
#   - neg_ll(tvec, param)
#   - func_T(tvec, param)
#   - K2operator(tvec, param, x, y)
#   - K2operatorAK2AT(tvec, param, A)
#   - K4operatorAABB(tvec, param, Q1, Q2)
#   - K3K3operatorAABBCC(tvec, param, Q1, Q2, Q3)
#   - K3K3operatorABCABC(tvec, param, Q1, Q2, Q3)
#   - K4operatorAABB_factored(tvec, param, A1, d1, A2, d2)
#   - K3K3operatorAABBCC_factored(tvec, param, A1, d1, A2, d2, A3, d3)
#   - K3K3operatorABCABC_factored(tvec, param, A1, d1, A2, d2, A3, d3)
#   - ineq_constraint(tvec, param)
#   - compute_analytic_tvec_hat(y, param)
# 
# Other objects:
#  - analytic_tvec_hat_func: function to compute tvec from y and parameters
#  - op_name: name of the operation (for call history)/for debugging: currently used in print method
#
# The 'initialize' method accepts these functions/objects or uses defaults.
# 
# Compulsory Methods:
# The class enforces the provision of five essential methods (K, K1, K2, K3operator, K4operator)
# by requiring their corresponding functions (K_func, K1_func, etc.) in the initialize method.
# 
# Optional Methods: 
# These methods like tilting_exponent, neg_ll, and func_T have default implementations.
# These defaults are used unless the user provides their own implementations by overriding the corresponding function.
# 
# EXTENSIBILITY:
#   'initialize' accepts '...' which allows adding new named methods
#   or overrides in future, without changing this class definition.
#   For instance, createCGF(...) can pass additional user-defined
#   methods that are stored or used in advanced ways.
# --------------------------------------------------------------------




#' @noRd
CGF <- R6::R6Class(
  classname = "CGF",
  
  # ----------------------------------------------------------
  # Private fields / methods
  # ----------------------------------------------------------
  private = list(
    
    # Core user-supplied method fields:
    K_func = NULL,
    K1_func = NULL,
    K2_func = NULL,
    K3operator_func = NULL,
    K4operator_func = NULL,
    
    
    # Optional method fields:
    ineq_constraint_func = NULL,
    ##### param_adaptor = NULL,    # We keep it here, but now set to identity by default
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
    
    # # ----------------------------------------------------------
    # # A helper method that used to adapt parameters,
    # # but now simply calls the passed-in func with param unmodified.
    # # We keep it so that the class structure remains unchanged.
    # # ----------------------------------------------------------
    # call_with_params = function(tvec, parameter_vector, func, ...) {
    #   # param_adaptor is now effectively identity
    #   p <- private$param_adaptor(parameter_vector)
    #   func(tvec, p, ...)
    # },
    
    
    # ----------------------------------------------------------
    # For optional methods, store either the user-supplied 
    # function or a default closure:
    # ----------------------------------------------------------
    
    # Tilting exponent term :  K(t) - sum(t_i * K1(t)_i)
    tilting_exponent = if (!is.null(private$tilting_exponent_func)) {
      private$tilting_exponent_func
    } else {
      function(tvec, parameter_vector)  {  # Default tilting_exponent = K - sum(tvec * K1)
        self$K(tvec, parameter_vector) - sum(tvec*self$K1(tvec, parameter_vector))
      }
    },
    
    
    
    
    
    neg_ll = if(!is.null(private$neg_ll_func)) {
      private$neg_ll_func(tvec, parameter_vector)
    } else {
      function(tvec, parameter_vector) {
        # Default neg_ll = 0.5*logdet(K2) + 0.5*log((2*pi)^m) - tilting_exponent
        K_val <- self$K(tvec, parameter_vector)
        K1_val <- self$K1(tvec, parameter_vector)
        te <- private$tilting_exponent(tvec, parameter_vector)
        K2_val <- self$K2(tvec, parameter_vector)
        val_logdet <- determinant(K2_val, logarithm = TRUE)$modulus
        res <- 0.5 * val_logdet + 0.5 * length(tvec)*log(2*pi) - te
        res
      }
    },
    
    
    #  The correction term of the first-order saddlepoint approximation to the log-likelihood
    func_T = if (!is.null(private$func_T_func)) {
      # Specify Q = K2^{-1} by computing an LDLT decomposition for K2
      # and using it to compute the A and d arguments to the _factored form of the operator methods
      # To save repeated inversion, this is performed once rather than delegated to the _factored_inverse methods
      # Note that with K2 = P^T L D L^T P, we have Q = K2^{-1} = P^{-1} L^{-1}^T D^{-1} L^{-1} P^{-1}^T
      # and since P is a permutation matrix P^T = P^{-1} and this reduces to
      # Q = (L^{-1} P)^T D^{-1} L^{-1} P, i.e., A=(L^{-1}P)^T = P^T (L^T)^{-1}
        private$func_T_func(tvec, parameter_vector)
      } else {
        function(tvec, parameter_vector) {
          K2_val <- self$K2(tvec, parameter_vector)
          K2_inv <- solve(K2_val)
          
          chol_K2_inv <- chol(K2_inv) # since K2_inv is symmetric positive definite
          diag_K2_inv <- diag(chol_K2_inv)
          d = diag_K2_inv^2
          # Create A with unit diagonals
          # A = L (from LDL) analog constructed from Cholesky
          A = t(chol_K2_inv) %*% diag(1/diag_K2_inv)
          
          # Evaluate required operators
          K4_AABB <- private$K4operatorAABB_factored(tvec, parameter_vector, A, d, A, d)
          K3K3_AABBCC <- private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A, d, A, d, A, d)
          K3K3_ABCABC <- private$K3K3operatorABCABC_factored(tvec, parameter_vector, A, d, A, d, A, d)
          K4_AABB/8 - K3K3_AABBCC/8 - K3K3_ABCABC/12
        }
      },
    
    
    
    
    #---------------------------------------------------------------------------
    # Alternative forms of the above operators where the matrices Qk (k=1,2 or k=1,2,3) are specified in factored form
    
    # Forms where the d-by-d symmetric matrix Qk is specified as Qk = Ak Dk Ak^T, where
    # Ak is d-by-r
    # Dk is r-by-r diagonal with diagonal entries given by the vector dk
    # This factorisation may arise as, for instance,
    # an LDLT factorisation, with Ak square and lower triangular (or permuted lower triangular)
    # an eigenvalue/eigenvector decomposition, with Ak the square orthogonal matrix whose columns are (right) eigenvectors of Qk
    # and dk the vector of eigenvalues of Qk
    # or in any other way.
    # In particular, d and r may be different, i.e., Ak need not be square (and the different Ak's may have different r values)
    # Note that this factorisation implies Qk(i,j) = sum_{m=1,...,r} dk(m)*Ak(i,m)*Ak(j,m)
    #---------------------------------------------------------------------------
    
    # K4operatorAABB factored version
    K4operatorAABB_factored =  if (!is.null(private$K4operatorAABB_factored_func)) {
      # Key identity:
      #   sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
      # = sum_{m1=1,...,r1} sum_{m2=1,...,r2} d1(m1)*d2(m2)*( sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)*A2(i4,m2) )
      # Complexity: r1*r2 calls to K4operator
        private$K4operatorAABB_factored_func(tvec, parameter_vector, A1, d1, A2, d2)
      } else {
          function(tvec, parameter_vector, A1, d1, A2, d2) {
            r1 <- length(d1)
            r2 <- length(d2)
            res <- 0
            for (m1 in seq_len(r1)) {
              for (m2 in seq_len(r2)) {
                res <- res + d1[m1]*d2[m2]*self$K4operator(tvec, parameter_vector, A1[,m1], A1[,m1], A2[,m2], A2[,m2])
              }
            }
            res
          }
      },
    

    K3K3operatorAABBCC_factored = function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
      # Key identity:
      #   sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
      # = sum_{m2=1,...,r2} d2(m2)*( sum_{m1=1,...,r1} d1(m1)*(sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)) )
      #                          *( sum_{m3=1,...,r3} d3(m3)*(sum_{i4,...i6=1,...,d} K3(i4,i5,i6)*A2(i4,m2)*A3(i5,m3)*A3(i6,m3)) )
      # Complexity: r2*(r1+r3) calls to K3operator
      if (!is.null(private$K3K3operatorAABBCC_factored_func)) {
        private$call_with_params(tvec, parameter_vector, private$K3K3operatorAABBCC_factored_func, A1, d1, A2, d2, A3, d3)
      } else {
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
    },
    
    
    
    K3K3operatorABCABC_factored = function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
      # Key identity:
      #   sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
      # = sum_{m1=1,...,r1} sum_{m2=1,...,r2} sum_{m3=1,...,r3} d1(m1)*d2(m2)*d3(m3)*
      #           ( sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A2(i2,m2)*A3(i3,m3) )^2
      # since the sum involving i4,i5,i6 is identical to the sum involving i1,i2,i3 except for the labelling of the variables of summation
      # Complexity: r1*r2*r3 calls to K3operator
      if (!is.null(private$K3K3operatorABCABC_factored_func)) {
        private$call_with_params(tvec, parameter_vector, private$K3K3operatorABCABC_factored_func, A1, d1, A2, d2, A3, d3)
      } else {
        r1 <- length(d1)
        r2 <- length(d2)
        r3 <- length(d3)
        
        message("The discrepancy option/compute.funcT has initiated a computation that may take a few moments...")
        
        
        # total_iterations <- r1 * r2 * r3
        # pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
        # iteration <- 0
        res <- 0
        for (m1 in seq_len(r1)) {
          for (m2 in seq_len(r2)) {
            for (m3 in seq_len(r3)) {
              val <- self$K3operator(tvec, parameter_vector, A1[,m1], A2[,m2], A3[,m3])
              res <- res + d1[m1]*d2[m2]*d3[m3]*(val*val)
              # iteration <- iteration + 1
              # setTxtProgressBar(pb, iteration)
            }
          }
        }
        # close(pb)
        res
      }
    },
    
    
    
    
    
    

    
    compute_analytic_tvec_hat_private = function(y, parameter_vector) {
      # no need to check has_analytic_tvec_hat(); this is done in the public method.
      stopifnot(is.numeric(y), !any(y <= 0))
      
      local_param <- private$param_adaptor(parameter_vector)
      private$analytic_tvec_hat_func(y, local_param)
    }
    
    
    
    
    
  ),
  
  
  # ----------------------------------------------------------
  # Public members
  # ----------------------------------------------------------
  public = list(
    
    call_history = NULL, # store an operation name or chain of calls
    
    # constructor
    initialize = function(K_func, K1_func, K2_func, K3operator_func, K4operator_func,
                          ineq_constraint_func = NULL,
                          param_adaptor = function(x) x,
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
      # Store user or default methods
      private$K_func <- K_func
      private$K1_func <- K1_func
      private$K2_func <- K2_func
      private$K3operator_func <- K3operator_func
      private$K4operator_func <- K4operator_func
      
      private$ineq_constraint_func <- ineq_constraint_func
      private$param_adaptor <- param_adaptor
      private$analytic_tvec_hat_func <- analytic_tvec_hat_func
      
      private$tilting_exponent_func <- tilting_exponent_func
      private$neg_ll_func <- neg_ll_func
      private$func_T_func <- func_T_func
      
      private$K4operatorAABB_func <- K4operatorAABB_func
      private$K3K3operatorAABBCC_func <- K3K3operatorAABBCC_func
      private$K3K3operatorABCABC_func <- K3K3operatorABCABC_func
      
      private$K4operatorAABB_factored_func <- K4operatorAABB_factored_func
      private$K3K3operatorAABBCC_factored_func <- K3K3operatorAABBCC_factored_func
      private$K3K3operatorABCABC_factored_func <- K3K3operatorABCABC_factored_func
      
      private$K2operator_func <- K2operator_func
      private$K2operatorAK2AT_func <- K2operatorAK2AT_func
      
      # Label or track the "operation" in call_history
      if (!is.character(op_name) ) stop("'operation' must of type character")
      
      if (!is.null(self$call_history)) {
        self$call_history <- c(self$call_history, op_name)
      } else {
        self$call_history <- op_name
      }
      
      
    },
    
    # Core mandatory CGF methods
    #' @description The CGF at point tvec.
    #' @param tvec Numeric vector at which to evaluate the CGF.
    #' @param parameter_vector Numeric vector of distribution parameters.
    # #' @return Numeric scalar value.
    K = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, private$K_func)
    },
    
    #' @description First derivative of the CGF.
    #' @param tvec Numeric vector.
    #' @param parameter_vector Numeric vector of parameters.
    # #' @return Numeric vector (the gradient of K at tvec).
    K1 = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, private$K1_func)
    },
    
    #' @description Second derivative (Hessian) of the CGF. 
    #' @param tvec Numeric vector.
    #' @param parameter_vector Numeric vector of parameters.
    # #' @return Numeric matrix.
    K2 = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, private$K2_func)
    },
    
    #' @description Third-order operator for the CGF.
    #' @param tvec Numeric vector.
    #' @param v1,v2,v3 Numeric vectors on which K3 acts.
    #' @param parameter_vector Numeric vector of parameters.
    # #' @return Numeric scalar representing the third-order operator application.
    K3operator = function(tvec, parameter_vector, v1, v2, v3) {
      private$call_with_params(tvec, parameter_vector, private$K3operator_func, v1, v2, v3)
    },
    
    #' @description Fourth-order operator for the CGF.
    #' @param tvec Numeric vector.
    #' @param v1,v2,v3,v4 Numeric vectors on which K4 acts.
    #' @param parameter_vector Numeric vector of parameters.
    # #' @return Numeric scalar representing the fourth-order operator application.
    K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) {
      private$call_with_params(tvec, parameter_vector, private$K4operator_func, v1, v2, v3, v4)
    },
    
    # -----------------------------------------------------------------------
    # Default implementations (available once K, K1, K2, K3operator, K4operator are provided)
    # -----------------------------------------------------------------------
    
    
    
    #' @description x^T K2 y operator
    #' @param tvec,x,y,parameter_vector As described.
    # #' @return Numeric scalar.
    K2operator = function(tvec, parameter_vector, x, y) {
      if (!is.null(private$K2operator_func)) {
        private$call_with_params(tvec, parameter_vector, private$K2operator_func, x, y)
      } else {
        # Default: x^T K2 y
        K2_val <- self$K2(tvec, parameter_vector)
        as.vector(t(x) %*% (K2_val %*% y))
      }
    },
    
    #' @description A K2 A^T operator
    #' @param tvec,A,parameter_vector As described.
    # #' @return Matrix
    K2operatorAK2AT = function(tvec, parameter_vector, A) {
      if (!is.null(private$K2operatorAK2AT_func)) {
        private$call_with_params(tvec, parameter_vector, private$K2operatorAK2AT_func, A)
      } else {
        # Default: A %*% K2 %*% t(A)
        K2_val <- self$K2(tvec, parameter_vector)
        A %*% K2_val %*% t(A)
      }
    },
    
    #' @description K4operatorAABB operator
    #' @param tvec,Q1,Q2,parameter_vector As described.
    # #' @return Numeric scalar.
    K4operatorAABB = function(tvec, parameter_vector, Q1, Q2) {
      if (!is.null(private$K4operatorAABB_func)) {
        private$call_with_params(tvec, parameter_vector, private$K4operatorAABB_func, Q1, Q2)
      } else {
        # K4operatorAABB: sum_i1,...i4 K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        private$K4operatorAABB_factored(tvec, parameter_vector, A1, d1, A1, d2) # Q1 == Q2
      }
    },
    
    #' @description K3K3operatorAABBCC operator
    #' @param tvec,Q1,Q2,Q3,parameter_vector As described.
    # #' @return Numeric scalar.
    K3K3operatorAABBCC = function(tvec, parameter_vector, Q1, Q2, Q3) {
      if (!is.null(private$K3K3operatorAABBCC_func)) {
        private$call_with_params(tvec, parameter_vector, private$K3K3operatorAABBCC_func, Q1, Q2, Q3)
      } else {
        # K3K3operatorAABBCC: sum over 6 indices with Q1,Q2,Q3
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1) # Q1 == Q2 == Q3
      }
    },
    
    #' @description K3K3operatorABCABC operator
    #' @param tvec,Q1,Q2,Q3,parameter_vector As described.
    # #' @return Numeric scalar.
    K3K3operatorABCABC = function(tvec, parameter_vector, Q1, Q2, Q3) {
      if(!is.null(private$K3K3operatorABCABC_func)) {
        private$call_with_params(tvec, parameter_vector, private$K3K3operatorABCABC_func, Q1, Q2, Q3)
      } else {
        # K3K3operatorABCABC: sum over 6 indices in a different pattern
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        private$K3K3operatorABCABC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1)
      }
    },
    
   
    
    #' @description Inequality constraints function.
    #' @param tvec,parameter_vector As described.
    # #' @return Numeric vector (default empty).
    ineq_constraint = function(tvec, parameter_vector) {
      if (!is.null(private$ineq_constraint_func)) {
        private$call_with_params(tvec, parameter_vector, private$ineq_constraint_func)
      } else {
        numeric(0)
      }
    },
    
    has_analytic_tvec_hat = function() !is.null(private$analytic_tvec_hat_func),
    analytic_tvec_hat = function(y, parameter_vector) {
      if (!self$has_analytic_tvec_hat()) return(NULL)
      private$compute_analytic_tvec_hat_private(y, parameter_vector)
    },
    
    
    print = function(...) {
      cat("<CGF Object>\n")
      if (!is.null(self$call_history)) {
        mapping_str <- paste(self$call_history, collapse = " -> ")
        cat("Used:", mapping_str, "\n")
      }
      cat("Class hierarchy:", paste(class(self), collapse = " -> "), "\n")
      invisible(self)
    },
    

    # Allows controlled access to certain private methods
    # Primarily for internal use
    .get_private_method = function(method_name) {
      if (!is.character(method_name) || length(method_name) != 1) stop("'method_name' must be a single character string.")
      if (!method_name %in% names(private)) stop(paste0("'", method_name, "' is not a private method in the CGF class."))
      # a whitelist of allowed private methods
      # add other allowed method names here
      allowed_methods <- c("neg_ll", "tilting_exponent", 
                           "K4operatorAABB_factored", 
                           "K3K3operatorAABBCC_factored", 
                           "K3K3operatorABCABC_factored",
                           "func_T"
                           )
      if (!(method_name %in% allowed_methods)) stop(paste0("Access to private method '", method_name, "' is not permitted."))
      
      method_ <- private[[method_name]]
      if (!is.function(method_)) stop(paste0("Private member '", method_name, "' is not a function."))
      method_
    },
    
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








#' Create a CGF object from user-defined functions
#'
#' 
#' @description
#' This creates an object of type CGF using user-supplied functions. You supply essential methods (`K`, `K1`, `K2`, etc.) plus any
#' optional overrides (e.g., `tilting_exponent` or `neg_ll`), and it returns
#' a `CGF` instance. This approach allows the user to define a CGF object using a set of functions
#' instead of requiring the user to define a subclass. This approach is more flexible
#' and avoids the need for inheritance in many cases.
#'
#'
#' @param K A function `K(tvec, parameter_vector) -> numeric scalar`.
#' @param K1 A function `K1(tvec, parameter_vector) -> numeric vector`.
#' @param K2 A function `K2(tvec, parameter_vector) -> numeric matrix`.
#' 
#' @param K3operator A function implementing the third-order operator.
#' @param K4operator A function implementing the fourth-order operator.
#' @param ineq_constraint Optional function for inequality constraints.
#' @param param_adaptor Optional function to adapt the parameter vector before
#'   calling the user methods. Defaults to the identity function.
#' 
#' @param analytic_tvec_hat_func Optional function for an analytic solution
#'   of the saddlepoint equation. If provided, call it via `cgf$analytic_tvec_hat(y, param)`.
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
#'
#' @return An object of class `CGF`.
#' @export
createCGF <- function(K, K1, K2, K3operator, K4operator, 
                      ineq_constraint = NULL,
                      param_adaptor = function(x) x,
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
  
  # Collect any additional methods passed via ...
  additional_methods <- list(...)
  
  # Merge user-supplied optional methods with additional methods
  # Additional methods take precedence in case of name conflicts
  all_optional_methods <- modifyList(user_optional_methods, additional_methods)
  
  
  do.call(CGF$new, c(
    list(
      K_func                = K,
      K1_func               = K1,
      K2_func               = K2,
      K3operator_func       = K3operator,
      K4operator_func       = K4operator,
      ineq_constraint_func  = ineq_constraint,
      param_adaptor    = param_adaptor,
      analytic_tvec_hat_func = analytic_tvec_hat_func,
      op_name = op_name
    ),
    all_optional_methods
  ))
  
}




