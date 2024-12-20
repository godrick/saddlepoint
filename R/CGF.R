# -------------------------------------------------------------------------
# This file defines an R class to represent a base CGF object.
#
#
# We provide a simplified R6 class "CGF" that:
# - Provides an interface: methods K, K1, K2, K3operator, K4operator, etc.
# - Implements defaults 
# - Methods that are distribution-specific (like K, K1, K2, etc.) are left
#   as placeholders to be implemented by subclasses.
# - "Operator" methods (like K2operator, K4operatorAABB_factored) rely on
#   K, K1, K2, K3operator, K4operator.
# -------------------------------------------------------------------------







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
#'
#' @export
CGF <- R6::R6Class(
  classname = "CGF",
  
  private = list(
    # A helper method that adapts parameters once and then calls the provided function
    call_with_params = function(tvec, parameter_vector, func, ...) {
      p <- self$param_adaptor(parameter_vector)
      func(tvec, p, ...)
    }
  ),
  
  public = list(
    K_func = NULL,
    K1_func = NULL,
    K2_func = NULL,
    K3operator_func = NULL,
    K4operator_func = NULL,
    ineq_constraint_func = NULL,
    param_adaptor = NULL,
    tilting_exponent_func = NULL,
    neg_ll_func = NULL,
    func_T_func = NULL,
    K2operator_func = NULL,
    K2operatorAK2AT_func = NULL,
    K4operatorAABB_func = NULL,
    K3K3operatorAABBCC_func = NULL,
    K3K3operatorABCABC_func = NULL,
    K4operatorAABB_factored_func = NULL,
    K3K3operatorAABBCC_factored_func = NULL,
    K3K3operatorABCABC_factored_func = NULL,
    
    
    #' @description Initialize the CGF object with user-defined functions and a parameter adaptor.
    #' @param K_func Function for K.
    #' @param K1_func Function for K1.
    #' @param K2_func Function for K2.
    #' @param K3operator_func Function for K3operator.
    #' @param K4operator_func Function for K4operator.
    #' @param ineq_constraint_func Function for inequality constraints.
    #' @param tilting_exponent_func Function for tilting exponent.
    #' @param neg_ll_func Function for negative log-likelihood.
    #' @param func_T_func Function for correction term.
    #' @param K2operator_func Function for K2operator.
    #' @param K2operatorAK2AT_func Function for K2operatorAK2AT.
    #' @param K4operatorAABB_func Function for K4operatorAABB.
    #' @param K3K3operatorAABBCC_func Function for K3K3operatorAABBCC.
    #' @param K3K3operatorABCABC_func Function for K3K3operatorABCABC.
    #' @param K4operatorAABB_factored_func Function for K4operatorAABB_factored.
    #' @param K3K3operatorAABBCC_factored_func Function for K3K3operatorAABBCC_factored.
    #' @param K3K3operatorABCABC_factored_func Function for K3K3operatorABCABC_factored.
    #' @param param_adaptor Function to adapt parameter_vector.
    initialize = function(K_func, K1_func, K2_func, K3operator_func, K4operator_func, 
                          ineq_constraint_func = NULL, 
                          param_adaptor = function(x) x,
                          tilting_exponent_func = NULL, 
                          neg_ll_func = NULL, 
                          func_T_func = NULL,
                          K2operator_func = NULL, 
                          K2operatorAK2AT_func = NULL, 
                          K4operatorAABB_func = NULL, 
                          K3K3operatorAABBCC_func = NULL, 
                          K3K3operatorABCABC_func = NULL,
                          K4operatorAABB_factored_func = NULL, 
                          K3K3operatorAABBCC_factored_func = NULL, 
                          K3K3operatorABCABC_factored_func = NULL
                          ) {
      # Assign functions
      self$K_func <- K_func
      self$K1_func <- K1_func
      self$K2_func <- K2_func
      self$K3operator_func <- K3operator_func
      self$K4operator_func <- K4operator_func
      self$ineq_constraint_func <- ineq_constraint_func
      self$tilting_exponent_func <- tilting_exponent_func
      self$neg_ll_func <- neg_ll_func
      self$func_T_func <- func_T_func
      self$K2operator_func <- K2operator_func
      self$K2operatorAK2AT_func <- K2operatorAK2AT_func
      self$K4operatorAABB_func <- K4operatorAABB_func
      self$K3K3operatorAABBCC_func <- K3K3operatorAABBCC_func
      self$K3K3operatorABCABC_func <- K3K3operatorABCABC_func
      self$K4operatorAABB_factored_func <- K4operatorAABB_factored_func
      self$K3K3operatorAABBCC_factored_func <- K3K3operatorAABBCC_factored_func
      self$K3K3operatorABCABC_factored_func <- K3K3operatorABCABC_factored_func
      self$param_adaptor <- param_adaptor
    },
    
    # -----------------------------------------------------------------------
    # Pure virtual-like methods that must be implemented by subclasses
    # The user implements specific distributions of interest (or parametric families of distributions)
    # as derived classes inheriting from the CGF class. 
    # -----------------------------------------------------------------------
    
    #' @description The CGF at point tvec.
    #' @param tvec Numeric vector at which to evaluate the CGF.
    #' @param parameter_vector Numeric vector of distribution parameters.
    #' @return Numeric scalar value.
    K = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, self$K_func)
    },
    
    #' @description First derivative of the CGF.
    #' @param tvec Numeric vector.
    #' @param parameter_vector Numeric vector of parameters.
    #' @return Numeric vector (the gradient of K at tvec).
    K1 = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, self$K1_func)
    },
    
    #' @description Second derivative (Hessian) of the CGF. 
    #' @param tvec Numeric vector.
    #' @param parameter_vector Numeric vector of parameters.
    #' @return Numeric matrix.
    K2 = function(tvec, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, self$K2_func)
    },
    
    #' @description Third-order operator for the CGF.
    #' @param tvec Numeric vector.
    #' @param v1,v2,v3 Numeric vectors on which K3 acts.
    #' @param parameter_vector Numeric vector of parameters.
    #' @return Numeric scalar representing the third-order operator application.
    K3operator = function(tvec, v1, v2, v3, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, self$K3operator_func, v1, v2, v3)
    },
    
    #' @description Fourth-order operator for the CGF.
    #' @param tvec Numeric vector.
    #' @param v1,v2,v3,v4 Numeric vectors on which K4 acts.
    #' @param parameter_vector Numeric vector of parameters.
    #' @return Numeric scalar representing the fourth-order operator application.
    K4operator = function(tvec, v1, v2, v3, v4, parameter_vector) {
      private$call_with_params(tvec, parameter_vector, self$K4operator_func, v1, v2, v3, v4)
    },
    
    # -----------------------------------------------------------------------
    # Default implementations (available once K, K1, K2, K3operator, K4operator are provided)
    # -----------------------------------------------------------------------
    
    
    #' @description Tilting exponent term : log numerator term in the saddlepoint approximation: K(t) - sum(t_i * K1(t)_i)
    #' @param tvec,parameter_vector As above.
    #' @return Numeric scalar.
    #' @export
    tilting_exponent = function(tvec, parameter_vector) {
      if (!is.null(self$tilting_exponent_func)) {
        private$call_with_params(tvec, parameter_vector, self$tilting_exponent_func)
      } else {
        # Default tilting_exponent = K - sum(tvec * K1)
        self$K(tvec, parameter_vector) - sum(tvec * self$K1(tvec, parameter_vector))
      }
    },
    
    
    #' @description The negative log-likelihood based on the saddlepoint approximation.
    #' @param tvec,parameter_vector As above.
    #' @return Numeric scalar.
    #' @export
    neg_ll = function(tvec, parameter_vector) {
      if(!is.null(self$neg_ll_func)) {
        private$call_with_params(tvec, parameter_vector, self$neg_ll_func)
      } else {
        # Default neg_ll = 0.5*logdet(K2) + 0.5*log((2*pi)^m) - tilting_exponent(
        K_val <- self$K(tvec, parameter_vector)
        K1_val <- self$K1(tvec, parameter_vector)
        te <- self$tilting_exponent(tvec, parameter_vector)
        K2_val <- self$K2(tvec, parameter_vector)
        val_logdet <- determinant(K2_val, logarithm = TRUE)$modulus
        res <- 0.5 * val_logdet + 0.5 * length(tvec)*log(2*pi) - te
        res
      }
    },
    

    #' @description The correction term of the first-order saddlepoint approximation to the log-likelihood
    #' @param tvec,parameter_vector As above.
    #' @return Numeric scalar.
    #' @export
    func_T = function(tvec, parameter_vector) {
      # Specify Q = K2^{-1} by computing an LDLT decomposition for K2
      # and using it to compute the A and d arguments to the _factored form of the operator methods
      # To save repeated inversion, this is performed once rather than delegated to the _factored_inverse methods
      # Note that with K2 = P^T L D L^T P, we have Q = K2^{-1} = P^{-1} L^{-1}^T D^{-1} L^{-1} P^{-1}^T
      # and since P is a permutation matrix P^T = P^{-1} and this reduces to
      # Q = (L^{-1} P)^T D^{-1} L^{-1} P, i.e., A=(L^{-1}P)^T = P^T (L^T)^{-1}
      if (!is.null(self$func_T_func)) {
        private$call_with_params(tvec, parameter_vector, self$func_T_func)
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
        K4_AABB <- self$K4operatorAABB_factored(tvec, A, d, A, d, parameter_vector)
        K3K3_AABBCC <- self$K3K3operatorAABBCC_factored(tvec, A, d, A, d, A, d, parameter_vector)
        K3K3_ABCABC <- self$K3K3operatorABCABC_factored(tvec, A, d, A, d, A, d, parameter_vector)
        
        K4_AABB/8 - K3K3_AABBCC/8 - K3K3_ABCABC/12
      }
    },
    
    #' @description x^T K2 y operator
    #' @param tvec,x,y,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K2operator = function(tvec, x, y, parameter_vector) {
      if (!is.null(self$K2operator_func)) {
        private$call_with_params(tvec, parameter_vector, self$K2operator_func, x, y)
      } else {
        # Default: x^T K2 y
        K2_val <- self$K2(tvec, parameter_vector)
        t(x) %*% (K2_val %*% y)
      }
    },
    
    #' @description A K2 A^T operator
    #' @param tvec,A,parameter_vector As described.
    #' @return Matrix
    #' @export
    K2operatorAK2AT = function(tvec, A, parameter_vector) {
      if (!is.null(self$K2operatorAK2AT_func)) {
        private$call_with_params(tvec, parameter_vector, self$K2operatorAK2AT_func, A)
      } else {
        # Default: A %*% K2 %*% t(A)
        K2_val <- self$K2(tvec, parameter_vector)
        A %*% K2_val %*% t(A)
      }
    },
    
    #' @description K4operatorAABB operator
    #' @param tvec,Q1,Q2,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K4operatorAABB = function(tvec, Q1, Q2, parameter_vector) {
      if (!is.null(self$K4operatorAABB_func)) {
        private$call_with_params(tvec, parameter_vector, self$K4operatorAABB_func, Q1, Q2)
      } else {
        # K4operatorAABB: sum_i1,...i4 K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        self$K4operatorAABB_factored(tvec, A1, d1, A1, d2, parameter_vector) # Q1 == Q2
      }
    },
    
    #' @description K3K3operatorAABBCC operator
    #' @param tvec,Q1,Q2,Q3,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K3K3operatorAABBCC = function(tvec, Q1, Q2, Q3, parameter_vector) {
      if (!is.null(self$K3K3operatorAABBCC_func)) {
        private$call_with_params(tvec, parameter_vector, self$K3K3operatorAABBCC_func, Q1, Q2, Q3)
      } else {
        # K3K3operatorAABBCC: sum over 6 indices with Q1,Q2,Q3
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        self$K3K3operatorAABBCC_factored(tvec, A1, d1, A1, d1, A1, d1, parameter_vector) # Q1 == Q2 == Q3
      }
    },
    
    #' @description K3K3operatorABCABC operator
    #' @param tvec,Q1,Q2,Q3,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K3K3operatorABCABC = function(tvec, Q1, Q2, Q3, parameter_vector) {
      if(!is.null(self$K3K3operatorABCABC_func)) {
        private$call_with_params(tvec, parameter_vector, self$K3K3operatorABCABC_func, Q1, Q2, Q3)
      } else {
        # K3K3operatorABCABC: sum over 6 indices in a different pattern
        chol_Q1 <- chol(Q1)
        diag_Q1 <- diag(chol_Q1)
        d1 <- diag_Q1^2
        A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
        self$K3K3operatorABCABC_factored(tvec, A1, d1, A1, d1, A1, d1, parameter_vector)
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
    #' @return Numeric scalar.
    #' @export
    K4operatorAABB_factored = function(tvec, A1, d1, A2, d2, parameter_vector) {
      # Key identity:
      #   sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
      # = sum_{m1=1,...,r1} sum_{m2=1,...,r2} d1(m1)*d2(m2)*( sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)*A2(i4,m2) )
      # Complexity: r1*r2 calls to K4operator
      if (!is.null(self$K4operatorAABB_factored_func)) {
        private$call_with_params(tvec, parameter_vector, self$K4operatorAABB_factored_func, A1, d1, A2, d2)
      } else {
        r1 <- length(d1)
        r2 <- length(d2)
        res <- 0
        for (m1 in seq_len(r1)) {
          for (m2 in seq_len(r2)) {
            res <- res + d1[m1]*d2[m2]*self$K4operator(tvec, A1[,m1], A1[,m1], A2[,m2], A2[,m2], parameter_vector)
          }
        }
        res
      }
    },
    
    #' @description K3K3operatorAABBCC factored version
    #' @param tvec,A1,d1,A2,d2,A3,d3,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K3K3operatorAABBCC_factored = function(tvec, A1, d1, A2, d2, A3, d3, parameter_vector) {
      # Key identity:
      #   sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
      # = sum_{m2=1,...,r2} d2(m2)*( sum_{m1=1,...,r1} d1(m1)*(sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)) )
      #                          *( sum_{m3=1,...,r3} d3(m3)*(sum_{i4,...i6=1,...,d} K3(i4,i5,i6)*A2(i4,m2)*A3(i5,m3)*A3(i6,m3)) )
      # Complexity: r2*(r1+r3) calls to K3operator
      if (!is.null(self$K3K3operatorAABBCC_factored_func)) {
        private$call_with_params(tvec, parameter_vector, self$K3K3operatorAABBCC_factored_func, A1, d1, A2, d2, A3, d3)
      } else {
        r1 <- length(d1)
        r2 <- length(d2)
        r3 <- length(d3)
        res <- 0
        for (m2 in seq_len(r2)) {
          factor1 <- 0
          for (m1 in seq_len(r1)) {
            factor1 <- factor1 + d1[m1]*self$K3operator(tvec, A1[,m1], A1[,m1], A2[,m2], parameter_vector)
          }
          factor2 <- 0
          for (m3 in seq_len(r3)) {
            factor2 <- factor2 + d3[m3]*self$K3operator(tvec, A2[,m2], A3[,m3], A3[,m3], parameter_vector)
          }
          res <- res + d2[m2]*factor1*factor2
        }
        res
      }
    },
    
    #' @description K3K3operatorABCABC factored version
    #' @param tvec,A1,d1,A2,d2,A3,d3,parameter_vector As described.
    #' @return Numeric scalar.
    #' @export
    K3K3operatorABCABC_factored = function(tvec, A1, d1, A2, d2, A3, d3, parameter_vector) {
      # Key identity:
      #   sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
      # = sum_{m1=1,...,r1} sum_{m2=1,...,r2} sum_{m3=1,...,r3} d1(m1)*d2(m2)*d3(m3)*
      #           ( sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A2(i2,m2)*A3(i3,m3) )^2
      # since the sum involving i4,i5,i6 is identical to the sum involving i1,i2,i3 except for the labelling of the variables of summation
      # Complexity: r1*r2*r3 calls to K3operator
      if (!is.null(self$K3K3operatorABCABC_func)) {
        private$call_with_params(tvec, parameter_vector, self$K3K3operatorABCABC_func, A1, d1, A2, d2, A3, d3)
      } else {
        r1 <- length(d1)
        r2 <- length(d2)
        r3 <- length(d3)
        res <- 0
        for (m1 in seq_len(r1)) {
          for (m2 in seq_len(r2)) {
            for (m3 in seq_len(r3)) {
              val <- self$K3operator(tvec, A1[,m1], A2[,m2], A3[,m3], parameter_vector)
              res <- res + d1[m1]*d2[m2]*d3[m3]*(val*val)
            }
          }
        }
        res
      }
    },
    
    #' @description Inequality constraints function.
    #' @param tvec,parameter_vector As described.
    #' @return Numeric vector (default empty).
    #' @export
    ineq_constraint = function(tvec, parameter_vector) {
      if (!is.null(self$ineq_constraint_func)) {
        private$call_with_params(tvec, parameter_vector, self$ineq_constraint_func)
      } else {
        numeric(0)
      }
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
#'
#' @return An object of class `CGF`.
#' @export
createCGF <- function(K, K1, K2, K3operator, K4operator, 
                      ineq_constraint = NULL,
                      param_adaptor = function(x) {x},
                      tilting_exponent = NULL,
                      neg_ll = NULL,
                      func_T = NULL,
                      K2operator = NULL,
                      K2operatorAK2AT = NULL,
                      K4operatorAABB = NULL,
                      K3K3operatorAABBCC = NULL,
                      K3K3operatorABCABC = NULL,
                      K4operatorAABB_factored = NULL,
                      K3K3operatorAABBCC_factored = NULL,
                      K3K3operatorABCABC_factored = NULL
                      ) {
  # Check required functions
  if (!is.function(K))  stop("K is required and must be a function")
  if (!is.function(K1)) stop("K1 is required and must be a function")
  if (!is.function(K2)) stop("K2 is required and must be a function")
  if (!is.function(K3operator)) stop("K3operator must be a function")
  if (!is.function(K4operator)) stop("K4operator must be a function")
  #### TODO: Add checks for other functions
  
  CGF$new(
    K_func = K, 
    K1_func = K1, 
    K2_func = K2, 
    K3operator_func = K3operator, 
    K4operator_func = K4operator, 
    ineq_constraint_func = ineq_constraint, 
    param_adaptor = param_adaptor,
    tilting_exponent_func = tilting_exponent, 
    neg_ll_func = neg_ll, 
    func_T_func = func_T,
    K2operator_func = K2operator, 
    K2operatorAK2AT_func = K2operatorAK2AT, 
    K4operatorAABB_func = K4operatorAABB, 
    K3K3operatorAABBCC_func = K3K3operatorAABBCC, 
    K3K3operatorABCABC_func = K3K3operatorABCABC,
    K4operatorAABB_factored_func = K4operatorAABB_factored_func, 
    K3K3operatorAABBCC_factored_func = K3K3operatorAABBCC_factored_func, 
    K3K3operatorABCABC_factored_func = K3K3operatorABCABC_factored_func
  )
}




