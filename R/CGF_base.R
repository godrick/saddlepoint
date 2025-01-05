# R/CGF_base.R
# Objects: CGF, createCGF

# ------------------------------------------------------------------------------
# Base CGF Class
#
# This class represents a cumulant generating function (CGF) and its related
# operators. Required methods (K, K1, K2, K3operator, K4operator) must be given.
# Optional methods (such as tilting_exponent, neg_ll, and various operator forms)
# can be supplied or will default to a basic implementation.
#
# The methods are:
# 1. K
# 2. K1
# 3. K2
# 4. K3operator
# 5. K4operator
# 6. ineq_constraint
# 7. param_adaptor
# 8. tilting_exponent
# 9. neg_ll
# 10. K4operatorAABB
# 11. K3K3operatorAABBCC
# 12. K3K3operatorABCABC
# 13. func_T
# 14. K4operatorAABB_factored
# 15. K3K3operatorAABBCC_factored
# 16. K3K3operatorABCABC_factored
# 17. K2operator
# 18. K2operatorAK2AT
# ------------------------------------------------------------------------------


# Compulsory Methods: 
# The class enforces the provision of five essential methods (K, K1, K2, K3operator, K4operator) 
# by requiring their corresponding functions (K_func, K1_func, etc.) in the initialize method.
#
# Optional Methods: Methods like tilting_exponent, neg_ll, and func_T have default implementations within the class. 
# These defaults are used unless the user provides their own implementations by overriding the 
# corresponding function.
#
# Additional Operator Methods: 
# Methods such as K4operatorAABB, K3K3operatorAABBCC, and K3K3operatorABCABC 
# also have default implementations, which can be overridden.
#
# Extensibility: 
# The initialize method accepts ..., 
# allowing users to pass additional named optional methods in the future 
# without modifying the class definition.





#' Base CGF Class
#'
#' @description
#' `CGF` provides a base class for implementing cumulant generating functions (CGFs) of various distributions.
#' 
#' Distributions of interest should define these methods:
#' - `K(tvec, parameter_vector)`: The CGF at point `tvec`.
#' - `K1(tvec, parameter_vector)`: The first tvec-derivative of the CGF.
#' - `K2(tvec, parameter_vector)`: The second-order tvec-gradient of the CGF.
#' - `K3operator(tvec, v1, v2, v3, parameter_vector)`: The third-order operator.
#' - `K4operator(tvec, v1, v2, v3, v4, parameter_vector)`: The fourth-order operator.
#'
#' This class also provides default implementations for methods like `tilting_exponent()`, `neg_ll()`, and `func_T()`.
#'
#' @section Methods:
#' - `initialize(K_func, K1_func, K2_func, K3operator_func, K4operator_func, parameters=list())`: Initialize a CGF object.
#' - `K()`, `K1()`, `K2()`, `K3operator()`, `K4operator()`: Must be provided.
#' - `tilting_exponent(tvec, parameter_vector)`: Computes the log numerator term in the saddlepoint approximation: K(t) - sum(t_i * K1(t)_i)
#' - `neg_ll(tvec, parameter_vector)`: Computes the saddlepoint negative log-likelihood.
#' - `func_T(tvec, parameter_vector)`: Computes the correction term of the first-order saddlepoint approximation to the log-likelihood
#' - `K2operator()`, `K2operatorAK2AT()`, `K4operatorAABB()`, `K3K3operatorAABBCC()`, `K3K3operatorABCABC()`: Provide various operator forms.
#' - `K4operatorAABB_factored()`, `K3K3operatorAABBCC_factored()`, `K3K3operatorABCABC_factored()`: Factored forms of operator methods.
#' - `ineq_constraint(tvec, parameter_vector)`: Returns inequality constraints (default empty).
#' @noRd
CGF <- R6::R6Class(
  classname = "CGF",
  
  private = list(
    K_func = NULL,
    K1_func = NULL,
    K2_func = NULL,
    K3operator_func = NULL,
    K4operator_func = NULL,
    
    ineq_constraint_func = NULL,
    param_adaptor = NULL,
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
    
    # A helper method that adapts parameters once and then calls the provided function
    call_with_params = function(tvec, parameter_vector, func, ...) {
      p <- private$param_adaptor(parameter_vector)
      func(tvec, p, ...)
    },
    
    
    #' @description Tilting exponent term : log numerator term in the saddlepoint approximation: K(t) - sum(t_i * K1(t)_i)
    #' @param tvec,parameter_vector As above.
    # #' @return Numeric scalar.
    tilting_exponent = function(tvec, parameter_vector) {
      if (!is.null(private$tilting_exponent_func)) {
        private$call_with_params(tvec, parameter_vector, private$tilting_exponent_func)
      } else {
        # Default tilting_exponent = K - sum(tvec * K1)
        self$K(tvec, parameter_vector) - sum(tvec * self$K1(tvec, parameter_vector))
      }
    },
    
    #' @description The negative log-likelihood based on the saddlepoint approximation.
    #' @param tvec,parameter_vector As above.
    # #' @return Numeric scalar.
    neg_ll = function(tvec, parameter_vector) {
      if(!is.null(private$neg_ll_func)) {
        private$call_with_params(tvec, parameter_vector, private$neg_ll_func)
      } else {
        # Default neg_ll = 0.5*logdet(K2) + 0.5*log((2*pi)^m) - tilting_exponent(
        K_val <- self$K(tvec, parameter_vector)
        K1_val <- self$K1(tvec, parameter_vector)
        te <- private$tilting_exponent(tvec, parameter_vector)
        K2_val <- self$K2(tvec, parameter_vector)
        val_logdet <- determinant(K2_val, logarithm = TRUE)$modulus
        res <- 0.5 * val_logdet + 0.5 * length(tvec)*log(2*pi) - te
        res
      }
    },
    
    
    #' @description The correction term of the first-order saddlepoint approximation to the log-likelihood
    #' @param tvec,parameter_vector As above.
    # #' @return Numeric scalar.
    func_T = function(tvec, parameter_vector) {
      # Specify Q = K2^{-1} by computing an LDLT decomposition for K2
      # and using it to compute the A and d arguments to the _factored form of the operator methods
      # To save repeated inversion, this is performed once rather than delegated to the _factored_inverse methods
      # Note that with K2 = P^T L D L^T P, we have Q = K2^{-1} = P^{-1} L^{-1}^T D^{-1} L^{-1} P^{-1}^T
      # and since P is a permutation matrix P^T = P^{-1} and this reduces to
      # Q = (L^{-1} P)^T D^{-1} L^{-1} P, i.e., A=(L^{-1}P)^T = P^T (L^T)^{-1}
      if (!is.null(private$func_T_func)) {
        private$call_with_params(tvec, parameter_vector, private$func_T_func)
      } else {
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
    
    #' @description K4operatorAABB factored version
    #' @param tvec,A1,d1,A2,d2,parameter_vector As described.
    # #' @return Numeric scalar.
    K4operatorAABB_factored = function(tvec, parameter_vector, A1, d1, A2, d2) {
      # Key identity:
      #   sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
      # = sum_{m1=1,...,r1} sum_{m2=1,...,r2} d1(m1)*d2(m2)*( sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)*A2(i4,m2) )
      # Complexity: r1*r2 calls to K4operator
      if (!is.null(private$K4operatorAABB_factored_func)) {
        private$call_with_params(tvec, parameter_vector, private$K4operatorAABB_factored_func, A1, d1, A2, d2)
      } else {
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
    
    #' @description K3K3operatorAABBCC factored version
    #' @param tvec,A1,d1,A2,d2,A3,d3,parameter_vector As described.
    # #' @return Numeric scalar.
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
    
    #' @description K3K3operatorABCABC factored version
    #' @param tvec,A1,d1,A2,d2,A3,d3,parameter_vector As described.
    #' @return Numeric scalar.
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
    
    
    
    
    
    #' @description Compute tvec using analytic_tvec_hat_func.
    compute_analytic_tvec_hat_private = function(y, parameter_vector) {
      if (!is.null(private$analytic_tvec_hat_func)) {
        stopifnot(is.numeric(y), !any(y <= 0))
        p <- private$param_adaptor(parameter_vector)
        private$analytic_tvec_hat_func(y, p)
      } else {
        # warning("analytic_tvec_hat is not available for this CGF object.")
        NULL
      }
    }
    
    
    
    
  ),
  
  public = list(
    
    # Initialize a CGF object
    # 
    # @param K_func,K1_func,K2_func,K3operator_func,K4operator_func Functions for required methods.
    # @param ineq_constraint_func Optional function for inequality constraints.
    # @param param_adaptor Optional function to adapt parameters (default identity).
    # @param tilting_exponent_func Optional tilting exponent function.
    # @param neg_ll_func Optional negative log-likelihood function.
    # @param func_T_func Optional correction term function.
    # @param K4operatorAABB_func,K3K3operatorAABBCC_func,K3K3operatorABCABC_func Optional operator functions.
    # @param K4operatorAABB_factored_func,K3K3operatorAABBCC_factored_func,K3K3operatorABCABC_factored_func Optional factored operator functions.
    # @param K2operator_func,K2operatorAK2AT_func Optional operator functions.
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
      # Assign functions
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
      
      
      if (!is.character(op_name) ) stop("'operation' must of type character")
      # Will store the history of calls or transformations.
      # # self$call_history <- op_name # 
      if (!is.null(self$call_history)) {
        self$call_history <- c(self$call_history, op_name)
      } else {
        self$call_history <- op_name
      }
      
      # # additional methods passed via ...
      # additional_methods <- list(...)
      # if (length(additional_methods) > 0) {
      #   for (name in names(additional_methods)) {
      #     self[[name]] <- additional_methods[[name]]
      #   }
      # }
    },
    
    # -----------------------------------------------------------------------
    # Pure virtual-like methods that must be implemented by subclasses
    # The user implements specific distributions of interest (or parametric families of distributions)
    # as derived classes inheriting from the CGF class. 
    # -----------------------------------------------------------------------
    
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
    
    
    analytic_tvec_hat = function(y, parameter_vector) {
      private$compute_analytic_tvec_hat_private(y, parameter_vector)
    },
    
    call_history = NULL, # Will store the history of calls or transformations.
    
    print = function(...) {
      cat("<CGF Object>\n")
      if (!is.null(self$call_history)) {
        mapping_str <- paste(self$call_history, collapse = " -> ")
        cat("Used:", mapping_str, "\n")
      }
      cat("Class hierarchy:", paste(class(self), collapse = " -> "), "\n")
      invisible(self)
    },
    

    # A public method to control access to specific private methods
    # Primarily for internal use, but can be used by the user if needed.
    # Hopefully, this discourages users from accessing these private methods.
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








#' A factory function to create a CGF object from user-defined functions
#'
#' @description
#' `createCGF()` creates an object that represents a CGF using user-supplied functions,
#' instead of requiring the user to define a subclass. This approach is more flexible
#' and avoids the need for inheritance in many cases.
#'
#' @param K A function implementing K(tvec, parameter_vector) returning a scalar.
#' @param K1 A function implementing K1(tvec, parameter_vector) returning a numeric vector.
#' @param K2 A function implementing K2(tvec, parameter_vector) returning a numeric matrix.
#' @param K3operator A function implementing ...
#' @param K4operator A function implementing ...
#' @param ineq_constraint (optional) If supplied, overrides the default `ineq_constraint`.
#' @param param_adaptor (optional) A function: param_adaptor(parameter_vector) returning adapted parameters.
#'                      Default is the identity function.
#' @param analytic_tvec_hat_func (Optional) function to compute tvec using analytical formula.
#' @param op_name A character string indicating the operation or transformation being performed.
#'                  Default is "UnnamedOperation".
#' @param tilting_exponent (optional) If supplied, overrides the default `tilting_exponent`.
#' @param neg_ll (optional) If supplied, overrides the default `neg_ll`.
#' @param func_T (optional) If supplied, overrides the default `func_T`.
#' @param K2operator (optional) If supplied, overrides the default `K2operator`.
#' @param K2operatorAK2AT (optional) If supplied, overrides the default `K2operatorAK2AT`.
#' @param K4operatorAABB (optional) If supplied, overrides the default `K4operatorAABB`.
#' @param K3K3operatorAABBCC (optional) If supplied, overrides the default `K3K3operator`.
#' @param K3K3operatorABCABC (optional) If supplied, overrides the default `K3K3operatorABCABC`.
#' @param K4operatorAABB_factored (optional) If supplied, overrides the default `K4operatorAABB_factored`.
#' @param K3K3operatorAABBCC_factored (optional) If supplied, overrides the default `K3K3operatorAABBCC_factored`.
#' @param K3K3operatorABCABC_factored (optional) If supplied, overrides the default `K3K3operatorABCABC_factored`.
#' @param ... Additional named optional methods that can be passed to the CGF object.
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




