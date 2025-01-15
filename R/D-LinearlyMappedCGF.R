# R/linearlyMappedCGF.R
# Object: linearlyMappedCGF
#### decide on lines - 77, 83
#
# A linearly mapped CGF for Y = A X, where X is the random vector of the CGF 'cgf'.
# 'matrix_A' can be either:
# (1) A fixed numeric matrix (dense or sparse).
# (2) A function param_vector -> numeric matrix (dense or sparse).
#
# If matrix_A is numeric but not sparse, it is converted to a sparse matrix.
# If matrix_A is already a function, it should itself handle any needed sparse conversion.
#
# If `iidReps = 1`, the code behaves exactly like a single linearly mapped CGF.
# If `iidReps > 1`, then we treat length(tvec) as `block_size * iidReps`,
# chunk it into `iidReps` blocks, and for each block we do:
#   block_t -> A^T block_t -> base_cgf's K, K1, etc.
# Then we sum (K, K3operator, K4operator) or create block-diagonal (K2) or
# piecewise concatenations (K1) as needed.
#


# # A helper to ensure we get a sparse matrix each time from A_fun
# get_sparse_A <- function(A_fun, param) {
#   A_current <- A_fun(param)
#   if (!inherits(A_current, "sparseMatrix")) {
#     A_current <- Matrix::Matrix(A_current, sparse = TRUE)
#   }
#   A_current
# }


#' @title CGF Object of a linearly mapped random variable (Y = A X)
#' 
#' @description
#' Creates a CGF object for the random vector \eqn{Y = A(\theta) \, X}, where
#' \eqn{X} is described by the input CGF `cgf`. The argument `matrix_A` can be:
#' 
#' - **A numeric matrix** (dense or sparse).
#' - **A function**: \eqn{\theta \mapsto A(\theta)} returning a numeric matrix.
#'
#' If `matrix_A` is a function, it is called for each invocation of the CGF
#' methods to retrieve the current matrix (allowing parameter-dependent transformations).
#' 
#' 
#' Internally, if \code{iidReps == 1} (the default), you get a single-block
#' linear mapping. If \code{iidReps > 1}, we build the single-block version
#' and then call \code{\link{iidReplicatesCGF}} to replicate that block logic.
#' 
#' When \code{iidReps > 1}, we treat the input dimension 
#' as \code{iidReps} repeated blocks. Each block is multiplied by the same matrix \eqn{A} 
#' and forwarded to the underlying `cgf`. 
#'
#'
#' @param cgf An object of class `CGF` for the base distribution \eqn{X}.
#' @param matrix_A Either a numeric matrix (dense or sparse), or a function:
#'   \code{function(param) -> numeric matrix}.
#' @param iidReps Integer. The number of IID replicate blocks to create. 
#'   Must be a positive integer. If \code{iidReps = 1}, it behaves as a single
#'   linearly mapped CGF with no repetition.
#' @param ... Additional arguments passed to \code{\link{createCGF}} or possibly
#'   to \code{\link{iidReplicatesCGF}} if \code{iidReps > 1}.
#'
#' @return A `CGF` object for \eqn{Y = A X}. 
#' @export
linearlyMappedCGF <- function(cgf, matrix_A, iidReps = 1, ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be an object inheriting from class 'CGF'.")
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1) stop("'iidReps' must be a positive integer.")
  
  #---------------------------------------------
  # # Convert matrix_A to a function param->sparse
  #---------------------------------------------
  is_matrix_A_function <- is.function(matrix_A)
  A_fun <- NULL
  if (!is_matrix_A_function) {
    # matrix_A is presumably a numeric matrix
    if (!is.matrix(matrix_A)) stop("'matrix_A' must be a numeric matrix or a function returning a matrix.")
                               if (!inherits(matrix_A, "sparseMatrix")) matrix_A <- Matrix::Matrix(matrix_A, sparse = TRUE)
    A_fun <- function(param) matrix_A 
  } else {
    # matrix_A is already a function
    A_fun <- function(param) {
      A_ <- matrix_A(param)
      if (!inherits(A_, "sparseMatrix")) A_ <- Matrix::Matrix(A_, sparse = TRUE)
      A_
    }
  }
  
  #---------------------------------------------
  # # Single-block linearlyMapped CGF
  # # Overrides for K, K1, K2, etc., where 'A_current' = get_sparse_A(parameter_vector)
  #---------------------------------------------
  
  # A helper to get a sparse matrix once
  ### this is not needed, but I'll keep it for now; might be useful for exta checks (central spot for extra logic)
  get_sparse_A <- function(param) A_fun(param)
  
  
  # Key identity: K_Y(t) = K_X(A^T t) (with t assumed to be a column vector)
  Kfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    cgf$K(t(A_current) %*% tvec, parameter_vector)
  }
  
  # Key identity: K_Y' = A K_X'
  K1fun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    A_current %*% cgf$K1(as.vector(t(A_current) %*% tvec), parameter_vector)
  }
  
  # Key identity: K_Y'' = A K_X'' A^T
  K2fun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    k2_base <- cgf$K2(as.vector(t(A_current) %*% tvec), parameter_vector)
    A_current %*% k2_base %*% t(A_current)
  }
  
  # Key identity: K_Y(t) - t^T K_Y'(t) = K_X(A^T t) - t^T A K_X'(A^T t) = K_X(A^T t) - (A^T t)^T K_X'(A^T t)
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  tiltingfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    tilting_exponent(as.vector(t(A_current) %*% tvec), parameter_vector)
  }
  
  # neg_ll: cgf's neg_ll will be used.
  # # neg_ll <- cgf$.get_private_method("neg_ll")
  # negllfun <- NULL
  
  # Key identity: x^T K_Y'' y = x^T A K_X'' A^T y = (A^T x)^T K_X'' A^T y
  K2operatorfun <- function(tvec, parameter_vector, x, y) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K2operator(as.vector(t(A_current) %*% tvec), 
                   parameter_vector,
                   as.vector(t(A_current) %*% x), 
                   as.vector(t(A_current) %*% y), 
                   )
  }
  
  # Returns B K_Y'' B^T as a function of the supplied (non-parameter) argument B
  # Key identity: B K_Y'' B^T = B A K_X'' A^T B^T = (B A) K_X'' (B A)^T
  K2operatorAK2ATfun <- function(tvec, parameter_vector, B) {
    A_current <- get_sparse_A(parameter_vector)
    B_A <- B %*% A_current
    cgf$K2operatorAK2AT(as.vector(t(A_current) %*% tvec), parameter_vector, B_A) %*% t(B_A)
  }
  
  K3operatorfun <- function(tvec, parameter_vector, v1, v2, v3) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K3operator(as.vector(t(A_current) %*% tvec), 
                   parameter_vector,
                   as.vector(t(A_current) %*% v1), 
                   as.vector(t(A_current) %*% v2), 
                   as.vector(t(A_current) %*% v3)
                   )
  }
  
  K4operatorfun <- function(tvec, parameter_vector, v1, v2, v3, v4) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K4operator(as.vector(t(A_current) %*% tvec), 
                   parameter_vector,
                   as.vector(t(A_current) %*% v1),
                   as.vector(t(A_current) %*% v2),
                   as.vector(t(A_current) %*% v3),
                   as.vector(t(A_current) %*% v4)
                   )
  }
  
  
  # All the operator forms involving matrices Q are equivalent to applying the same method for BaseCGF with Q_inner = A^T Q A
  K4operatorAABBfun <- function(tvec, parameter_vector, Q1, Q2) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    cgf$K4operatorAABB(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner)
  }
  
  K3K3operatorAABBCCfun <- function(tvec, parameter_vector, Q1, Q2, Q3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    # Q3_inner <- tA %*% Q3 %*% A_current
    cgf$K3K3operatorAABBCC(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner, Q1_inner)
  }
  
  K3K3operatorABCABCfun <- function(tvec, parameter_vector, Q1, Q2, Q3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    # Q3_inner <- tA %*% Q3 %*% A_current
    cgf$K3K3operatorABCABC(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner, Q1_inner)
  }
  
  #### We avoid the factored forms for now (avoiding the potentially expensive loops)
  func_Tfun <- function(tvec, parameter_vector) {
                         Q <- solve(K2fun(tvec, parameter_vector))
    K3K3operatorABCABC_val <- K3K3operatorABCABCfun(tvec, parameter_vector, Q, Q, Q)
    K3K3operatorAABBCC_val <- K3K3operatorAABBCCfun(tvec, parameter_vector, Q, Q, Q)
        K4operatorAABB_val <- K4operatorAABBfun(tvec, parameter_vector, Q, Q)
    K4operatorAABB_val/8 - K3K3operatorAABBCC_val/8 - K3K3operatorABCABC_val/12
  }
  
  
  # For the factored forms where Q = B D B^T and D has diagonal vector d, note that Q_inner = A^T Q A = (A^T B) D (A^T B)^T
  # Note about sizes: if A is n-by-m then B is n-by-r for some r, and A^T B is m-by-r
  K4operatorAABB_factored <- cgf$.get_private_method("K4operatorAABB_factored")
  K4operatorAABB_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    K4operatorAABB_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2)
  }
  
  K3K3operatorAABBCC_factored <- cgf$.get_private_method("K3K3operatorAABBCC_factored")
  K3K3operatorAABBCC_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2, B3, d3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    B3_inner <- tA %*% B3
    K3K3operatorAABBCC_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2, B3_inner, d3)
  }
  
  K3K3operatorABCABC_factored <- cgf$.get_private_method("K3K3operatorABCABC_factored")
  K3K3operatorABCABC_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2, B3, d3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    B3_inner <- tA %*% B3
    K3K3operatorABCABC_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2, B3_inner, d3)
  }
  
  # inequality constraints for the transformed variable Y = A * X are the same as those
  # for the original variable X, evaluated at the transformed input A.transpose() * tvec.
  ineq_constraintfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$ineq_constraint(as.vector(t(A_current) %*% tvec), parameter_vector)
  }
  
  
  # ------------------------------------------------------------------
  # # Build the new mapped CGF using createCGF
  # ------------------------------------------------------------------
  mapped_cgf <- createCGF(
    K = Kfun, 
    K1 = K1fun, 
    K2 = K2fun, 
    K3operator = K3operatorfun, 
    K4operator = K4operatorfun, 
    ineq_constraint = ineq_constraintfun,
    param_adaptor = function(x) x,
    analytic_tvec_hat_func = NULL,
    tilting_exponent = tiltingfun,
    # neg_ll = negllfun,
    func_T = func_Tfun,
    K4operatorAABB = K4operatorAABBfun,
    K3K3operatorAABBCC = K3K3operatorAABBCCfun,
    K3K3operatorABCABC = K3K3operatorABCABCfun,
    K4operatorAABB_factored = K4operatorAABB_factoredfun,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factoredfun,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factoredfun,
    K2operator = K2operatorfun,
    K2operatorAK2AT = K2operatorAK2ATfun,
    op_name = c(cgf$call_history, "linearlyMappedCGF"),
    ...
  )
  
  #---------------------------------------------
  # # If iidReps == 1 => done. Otherwise wrap
  # # with iidReplicatesCGF().
  #---------------------------------------------
  if (iidReps == 1) return(mapped_cgf)
  
  # For iidReps > 1, we call an external aggregator that replicates blocks.
  # Note: pass the single-block cgf to 'iidReplicatesCGF', which will 
  # handle the chunking logic for all methods.
  
  iidReplicatesCGF(mapped_cgf, iidReps = iidReps, ...)
}
