# R/IIDReplicatesCGF.R
# Objects: iidReplicatesCGF
# Areas to review are marked with `####`
##### For now we do not unify the two cases (iidReps and block_size)
##### So any change in one case should be reflected in the other case as well


# #' For each operator (\eqn{K, K1, K2}, etc.), \code{iidReplicatesCGF()}:
# #' \itemize{
# #'   \item Splits \code{tvec} (and any corresponding vectors \code{x,y,v1,...}) into sub-blocks.
# #'   \item Calls the original CGF operator on each sub-block.
# #'   \item Combines results by summation (for \eqn{K}, \eqn{K3}, \eqn{K4}), or by
# #'         forming block-diagonal (for \eqn{K2} and similar).
# #' }
# #'
# #' Effectively, if the original CGF is compatible with dimension \eqn{d}, and \code{iidReps=m},
# #' the resulting CGF will be compatible with dimension \eqn{m d}, chunked into \eqn{m} blocks each
# #' of size \eqn{d}.
# #'


#------------------------------------------------------------
# Helper functions 
#------------------------------------------------------------
# a helper for chunking
chunkIndices <- function(i, block_size) {
  seq.int((i - 1)*block_size + 1, i*block_size)
}


#' @title Replicate a CGF over multiple i.i.d blocks or with a fixed block size
#'
#' @description
#' Extends a given CGF to handle multiple IID blocks based on either `iidReps` 
#' or a fixed `block_size`. Exactly one of these parameters must be provided.
#'
#' @param cgf A CGF object. 
#' @param iidReps Either NULL or a positive integer specifying the number of IID blocks.
#' @param block_size Either NULL or a positive integer specifying the size of each block.
#' @param ... Additional arguments passed to \code{\link{createCGF}} in case
#'   you want to override any methods.
#'
#' @return A new CGF object that operates on length-\eqn{m d} input vectors, chunked
#'   into \eqn{m} blocks, each block of length \eqn{d}.
#'
#' @details
#' - If `iidReps` is provided, the input vector is split into `iidReps` equal-sized blocks.
#' - If `block_size` is provided, the number of blocks is determined by the length of the input vector.
#' - Exactly one of `iidReps` or `block_size` should be non-NULL.
#'
#' @examples
#' # Suppose cgf is dimension=3, and we want 5 blocks => dimension=15
#' #   aggregator <- iidReplicatesCGF(cgf, iidReps=5)
#'
#' @export
iidReplicatesCGF <- function(cgf, iidReps = NULL, block_size = NULL, ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object.")
  if (is.null(iidReps) && is.null(block_size)) stop("Specify either iidReps or block_size.")
  if (!is.null(iidReps) && !is.null(block_size)) stop("Specify only one of iidReps or block_size.")
  
  if (!is.null(iidReps)) {
    if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 2) {
      stop("'iidReps' must be an integer >= 2.")
    }
  }
  
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1 || block_size < 2) {
      stop("'block_size' must be an integer >= 2.")
    }
  }
  
  
  # fetch some private methods from the base CGF
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  neg_ll <- cgf$.get_private_method("neg_ll")
  func_T <- cgf$.get_private_method("func_T")
  
  # ------------------------------------------------------------------
  # First we handle the iidReps case
  # ------------------------------------------------------------------
  
  if(!is.null(iidReps)){
    # K => sum over blocks
    Kfun <- function(tvec, param) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K(tvec[idx], param)
      }
      total
    }
    
    # ------------------------------------------------------------------
    # #### NOTE:
    # The line `out_ <- numeric(n)` currently triggers an error during
    # Tape recording. To address this issue, there are two potential solutions:
    #
    # 1. Utilize the ADoverload package to overload the `[<-` operator.
    # 2. Modify the line to create an explicit dependency on a parameter,
    #    such as `out_ <- numeric(n) * param[1]`.
    #
    # For the time being, we use the second approach to avoid any
    # complications with operator overloading.
    # ------------------------------------------------------------------
    # K1 => piecewise concatenation
    K1fun <- function(tvec, param) {
      # "[<-" <- RTMB::ADoverload("[<-")
      n <- length(tvec)
      block_size <- n / iidReps
      out_ <- numeric(n) * param[1]  #### Modified to depend on `param`
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        out_[idx] <- cgf$K1(tvec[idx], param)
      }
      out_
    }
    
    # K2 => block-diagonal
    K2fun <- function(tvec, param) {
      n <- length(tvec)
      if (n %% iidReps != 0) stop("The length of 'tvec' must be divisible by 'iidReps'.")
      
      block_size <- n / iidReps
      big_mat <- Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        big_mat[idx, idx] <- cgf$K2(tvec[idx], param)
      }
      big_mat
    }
    
    # tilting_exponent => sum
    tiltingfun <- function(tvec, param) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + tilting_exponent(tvec[idx], param)
      }
      total
    }
    
    # neg_ll => sum
    negllfun <- function(tvec, param) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + neg_ll(tvec[idx], param)
      }
      total
    }
    
    # K2operator => sum
    K2operatorfun <- function(tvec, param, x, y) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K2operator(tvec[idx], param, x[idx], y[idx])
      }
      total
    }
    
    # K2operatorAK2AT => block-diagonal
    K2operatorAK2ATfun <- function(tvec, param, B) {
      n <- length(tvec)
      block_size <- n / iidReps
      big_mat <- Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
      for (i in 1:iidReps) {
        idx  <- chunkIndices(i, block_size)
        subB <- B[idx, idx, drop=FALSE]
        big_mat[idx, idx] <- cgf$K2operatorAK2AT(tvec[idx], param, subB)
      }
      big_mat
    }
    
    # K3operator => sum
    K3operatorfun <- function(tvec, param, v1, v2, v3) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
      }
      total
    }
    
    # K4operator => sum
    K4operatorfun <- function(tvec, param, v1, v2, v3, v4) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
      }
      total
    }
    
    # K4operatorAABB => sum
    K4operatorAABBfun <- function(tvec, param, Q1, Q2) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        # slice out sub-block of Q1, Q2
        ### (Verify?) Q1 is block-diagonal of size (n·iidReps) × (n·iidReps). The relevant sub-block is n×n.
        Q1sub <- Q1[idx, idx, drop=FALSE]
        # Q2sub <- Q2[idx, idx, drop=FALSE]
        total <- total + cgf$K4operatorAABB(tvec[idx], param, Q1sub, Q1sub)
      }
      total
    }
    
    # K3K3operatorAABBCC => sum
    K3K3operatorAABBCCfun <- function(tvec, param, Q1, Q2, Q3) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        Q1sub <- Q1[idx, idx, drop=FALSE]
        # Q2sub <- Q2[idx, idx, drop=FALSE]
        # Q3sub <- Q3[idx, idx, drop=FALSE]
        total <- total + cgf$K3K3operatorAABBCC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
      }
      total
    }
    
    # K3K3operatorABCABC => sum
    K3K3operatorABCABCfun <- function(tvec, param, Q1, Q2, Q3) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        Q1sub <- Q1[idx, idx, drop=FALSE]
        # Q2sub <- Q2[idx, idx, drop=FALSE]
        # Q3sub <- Q3[idx, idx, drop=FALSE]
        total <- total + cgf$K3K3operatorABCABC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
      }
      total
    }
    
    # func_T => sum
    func_Tfun <- function(tvec, param) {
      block_size <- length(tvec) / iidReps
      total <- 0
      for (i in 1:iidReps) {
        idx <- chunkIndices(i, block_size)
        total <- total + func_T(tvec[idx], param)
      }
      total
    }
    
    # ineq_constraint => concatenation
    ineq_constraintfun <- function(tvec, param) {
      block_size <- length(tvec) / iidReps
      
      # call once for length
      idx <- chunkIndices(1, block_size)
      first_val <- cgf$ineq_constraint(tvec[idx], param)
      L <- length(first_val)
      
      result <- numeric(L * iidReps)*param[1]  #### Modified to depend on `param`
      result[1:L] <- first_val
      for (i in 2:iidReps) {
        idx <- chunkIndices(i, block_size)
        val <- cgf$ineq_constraint(tvec[idx], param)
        # if (length(val) != L) stop("ineq_constraint changed length across blocks?")
        start_ <- (i-1)*L + 1
        result[start_:(i*L)] <- val
      }
      result
    }
    
    
    # for the analytic_tvec_hat_func:
    # We'll do a chunk approach if cgf$analytic_tvec_hat() is non-NULL:
    # e.g. chunk y => pass each chunk to cgf$analytic_tvec_hat => combine?
    #### Check if this doesn't make sense, (default to NULL if that's the case)
    if (cgf$has_analytic_tvec_hat() ) {
      # define a chunk aggregator
      analytic_tvec_hat_func <- function(y, param) {
        block_size <- length(y) / iidReps
        # we chunk y and call base_analytic_hat
        out_ <- numeric(length(y))*param[1]  #### Modified to depend on `param`
        for (i in 1:iidReps) {
          idx <- chunkIndices(i, block_size)
          out_[idx] <- cgf$analytic_tvec_hat(y[idx], param)
        }
        out_
      }
    } else {
      # If the base CGF had no valid function, just return NULL
      analytic_tvec_hat_func <- NULL
    }
  } else{
    # ------------------------------------------------------------------
    # Next we handle the block_size case
    # ------------------------------------------------------------------
    Kfun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K(tvec[idx], param)
      }
      total
    }
    
    K1fun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      n <- length(tvec)
      nBlocks <- n / block_size
      out_ <- numeric(n) * param[1]  #### Modified to depend on `param`
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        out_[idx] <- cgf$K1(tvec[idx], param)
      }
      out_
    }
    
    # K2 => block-diagonal
    K2fun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      n <- length(tvec)
      nBlocks <- n / block_size
      big_mat <- Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        big_mat[idx, idx] <- cgf$K2(tvec[idx], param)
      }
      big_mat
    }
    
    # tilting_exponent => sum
    tiltingfun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + tilting_exponent(tvec[idx], param)
      }
      total
    }
    
    # neg_ll => sum
    negllfun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + neg_ll(tvec[idx], param)
      }
      total
    }
    
    # K2operator => sum
    K2operatorfun <- function(tvec, param, x, y) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K2operator(tvec[idx], param, x[idx], y[idx])
      }
      total
    }
    
    # K2operatorAK2AT => block-diagonal
    K2operatorAK2ATfun <- function(tvec, param, B) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      n <- length(tvec)
      nBlocks <- n / block_size
      big_mat <- Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
      for (i in 1:nBlocks) {
        idx  <- chunkIndices(i, block_size)
        subB <- B[idx, idx, drop=FALSE]
        big_mat[idx, idx] <- cgf$K2operatorAK2AT(tvec[idx], param, subB)
      }
      big_mat
    }
    
    # K3operator => sum
    K3operatorfun <- function(tvec, param, v1, v2, v3) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
      }
      total
    }
    
    # K4operator => sum
    K4operatorfun <- function(tvec, param, v1, v2, v3, v4) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + cgf$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
      }
      total
    }
    
    # K4operatorAABB => sum
    K4operatorAABBfun <- function(tvec, param, Q1, Q2) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        # slice out sub-block of Q1, Q2
        ### (Verify?) Q1 is block-diagonal of size (n·iidReps) × (n·iidReps). The relevant sub-block is n×n.
        Q1sub <- Q1[idx, idx, drop=FALSE]
        # Q2sub <- Q2[idx, idx, drop=FALSE]
        total <- total + cgf$K4operatorAABB(tvec[idx], param, Q1sub, Q1sub)
      }
      total
    }
    
    # K3K3operatorAABBCC => sum
    K3K3operatorAABBCCfun <- function(tvec, param, Q1, Q2, Q3) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        Q1sub <- Q1[idx, idx, drop=FALSE]
        total <- total + cgf$K3K3operatorAABBCC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
      }
      total
    }
    
    # K3K3operatorABCABC => sum
    K3K3operatorABCABCfun <- function(tvec, param, Q1, Q2, Q3) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        Q1sub <- Q1[idx, idx, drop=FALSE]
        # Q2sub <- Q2[idx, idx, drop=FALSE]
        # Q3sub <- Q3[idx, idx, drop=FALSE]
        total <- total + cgf$K3K3operatorABCABC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
      }
      total
    }
    
    # func_T => sum
    func_Tfun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      total <- 0
      for (i in 1:nBlocks) {
        idx <- chunkIndices(i, block_size)
        total <- total + func_T(tvec[idx], param)
      }
      total
    }
    
    # ineq_constraint => concatenation
    ineq_constraintfun <- function(tvec, param) {
      if (length(tvec) %% block_size != 0) stop("Length of tvec must be divisible by block_size.")
      nBlocks <- length(tvec) / block_size
      
      # call once for length
      idx <- chunkIndices(1, block_size)
      first_val <- cgf$ineq_constraint(tvec[idx], param)
      L <- length(first_val)
      
      result <- numeric(L * nBlocks)*param[1]  #### Modified to depend on `param`
      result[1:L] <- first_val
      for (i in 2:nBlocks) {
        idx <- chunkIndices(i, block_size)
        val <- cgf$ineq_constraint(tvec[idx], param)
        start_ <- (i-1)*L + 1
        result[start_:(i*L)] <- val
      }
      result
    }
    
    
    if ( cgf$has_analytic_tvec_hat()  ) {
      # define a chunk aggregator
      analytic_tvec_hat_func <- function(y, param) {
        if (length(y) %% block_size != 0) stop("Length of y must be divisible by block_size.")
        nBlocks <- length(y) / block_size

        # we chunk y and call base_analytic_hat
        out_ <- numeric(length(y))*param[1]  #### Modified to depend on `param`
        for (i in 1:nBlocks) {
          idx <- chunkIndices(i, block_size)
          out_[idx] <- cgf$analytic_tvec_hat(y[idx], param)
        }
        out_
      }
    } else {
      # If the base CGF had no valid function, just return NULL
      analytic_tvec_hat_func <- NULL
    }
    
  }
  
  
  # Build a new CGF object
  replicated_cgf <- createCGF(
          K = Kfun, 
          K1 = K1fun, 
          K2 = K2fun, 
          K3operator = K3operatorfun, 
          K4operator = K4operatorfun, 
          ineq_constraint = ineq_constraintfun,
          param_adaptor = function(x) x,
          analytic_tvec_hat_func = analytic_tvec_hat_func,
          tilting_exponent = tiltingfun,
          neg_ll = negllfun,
          func_T = func_Tfun,
          K4operatorAABB = K4operatorAABBfun,
          K3K3operatorAABBCC = K3K3operatorAABBCCfun,
          K3K3operatorABCABC = K3K3operatorABCABCfun,
          # K4operatorAABB_factored = K4operatorAABB_factoredfun,
          # K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factoredfun,
          # K3K3operatorABCABC_factored = K3K3operatorABCABC_factoredfun,
          K2operator = K2operatorfun,
          K2operatorAK2AT = K2operatorAK2ATfun,
          op_name = c(cgf$call_history, "iidReplicatesCGF"),
          ...
  )
  
  replicated_cgf
}






