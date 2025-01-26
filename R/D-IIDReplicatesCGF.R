# R/IIDReplicatesCGF.R
# Objects: iidReplicatesCGF
# Areas to review are marked with `####`
##### For now we do not unify the two cases (iidReps and block_size)
##### So any change in one case should be reflected in the other case as well



###############################################################################
# NOTE:
#
# In this code, we refer to functions like K(), K1(), K2(), etc. as "aggregators."
# Each of these functions takes an input vector `tvec` and combines (sums or
# concatenates) results across multiple blocks (i.i.d. or fixed size). In other
# words, these functions "aggregate" the block-wise results from the underlying
# base CGF into a final output. That's why we call them "aggregators."
#
# Based on the arguments 'iidReps' and 'block_size', we may have three scenarios:
#
# (A) A single aggregator function with minimal branching:
#     - Store 'iidReps', 'block_size', etc. in closure variables.
#     - Do a small if/else for each aggregator (K, K1, etc.) to detect scenario.
#     - Pros: Little code duplication.
#     - Cons: Some branching in every call, and the code can get a bit messy.
#
# (B) Completely separate aggregator definitions for each scenario:
#     - One set of aggregator functions if (iidReps && block_size),
#       one set if only iidReps, another if only block_size.
#     - Pros: Zero branching inside aggregator calls.
#     - Cons: Code duplication across scenarios (harder to maintain).
#
# (C) Helper functions (get_nBlocks, get_blockSize) that do scenario checks:
#     - Each aggregator just calls these helpers, then loops over blocks.
#     - Pros: Easier to maintain than (B); simpler than (A) because the scenario logic is centralized
#       in the helpers.
#     - Cons: Slight overhead of helper function calls. Not sure if this can be a performance issue.
#
# ==> For now, we chose Option C
###############################################################################





#------------------------------------------------------------
# Helper functions 
#------------------------------------------------------------
# a helper for chunking
chunkIndices <- function(i, block_size) {
  seq.int((i - 1)*block_size + 1, i*block_size)
}





# ------------------------------------------------------------------
# Internal function; the actual implementation
#   iidReps and block_size will have been by the calling function.
# ------------------------------------------------------------------
.iidReplicatesCGF_internal <- function(cgf, iidReps, block_size) {
  useBoth   <- (!is.null(iidReps) && !is.null(block_size))
  onlyIID   <- (!is.null(iidReps) &&  is.null(block_size))
  onlyBlock <- ( is.null(iidReps) && !is.null(block_size))
  
  # ------------------------------------------------------------------
  # Helper functions that figure out nBlocks and blockSize on the fly
  #   depending on whether iidReps, block_size, or both are supplied.
  # ------------------------------------------------------------------
  get_nBlocks <- function(n) { # n here is the length of tvec
    if (onlyIID) {
      if (n %% iidReps != 0) {
        stop("Length of tvec must be divisible by 'iidReps'.")
      }
      return(iidReps) # number of blocks
    } 
    if (useBoth) {
      # Must match exactly block_size * iidReps
      if (n != block_size * iidReps) {
        stop("When both 'iidReps' and 'block_size' are given, ",
             "length(tvec) must be block_size * iidReps.")
      }
      return(iidReps)   
    } else {
      # onlyBlock
      if (n %% block_size != 0) {
        stop("Length of tvec must be divisible by 'block_size'.")
      }
      return(n / block_size)
    }
  }
  
  get_blockSize <- function(n) {
    if (onlyIID) {
      if (n %% iidReps != 0) {
        stop("Length of tvec must be divisible by 'iidReps'.")
      }
      # block size is n / iidReps
      return(n / iidReps)
    }
    if (useBoth) {
      # Must match exactly block_size * iidReps
      if (n != block_size * iidReps) {
        stop("When both 'iidReps' and 'block_size' are given, ",
             "length(tvec) must be block_size * iidReps.")
      }
      return(block_size)
    } else {
      # onlyBlock
      if (n %% block_size != 0) {
        stop("Length of tvec must be divisible by 'block_size'.")
      }
      return(block_size)
    }
  }
  
  
  # fetch some private methods from the base CGF
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  neg_ll <- cgf$.get_private_method("neg_ll")
  func_T <- cgf$.get_private_method("func_T")
  
  # ------------------------------------------------------------------
  # Now all methods in a unified manner
  # ------------------------------------------------------------------
  
  # K => sum over blocks
  Kfun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)  # how many blocks?
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
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
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    out_ <- numeric(n) * param[1]  # tie to 'param' to avoid tape error
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      out_[idx] <- cgf$K1(tvec[idx], param)
    }
    out_
  }
  
  # K2 => block-diagonal
  K2fun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    big_mat <- Matrix::Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      big_mat[idx, idx] <- cgf$K2(tvec[idx], param)
    }
    big_mat
  }
  
  # tilting_exponent => sum
  tiltingfun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + tilting_exponent(tvec[idx], param)
    }
    total
  }
  
  # neg_ll => sum
  negllfun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + neg_ll(tvec[idx], param)
    }
    total
  }
  
  
  # func_T => sum
  func_Tfun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + func_T(tvec[idx], param)
    }
    total
  }
  
  # K2operator => sum
  K2operatorfun <- function(tvec, param, x, y) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + cgf$K2operator(tvec[idx], param, x[idx], y[idx])
    }
    total
  }
  
  
  # K2operatorAK2AT => block-diagonal
  K2operatorAK2ATfun <- function(tvec, param, B) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    big_mat <- Matrix::Matrix(0, nrow = n, ncol = n) * param[1] #### decide on sparse?
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      subB <- B[idx, idx, drop = FALSE]
      big_mat[idx, idx] <- cgf$K2operatorAK2AT(tvec[idx], param, subB)
    }
    big_mat
  }
  
  # K3operator => sum
  K3operatorfun <- function(tvec, param, v1, v2, v3) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + cgf$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
    }
    total
  }
  
  # K4operator => sum
  K4operatorfun <- function(tvec, param, v1, v2, v3, v4) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      total <- total + cgf$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
    }
    total
  }

  
  # K4operatorAABB => sum
  K4operatorAABBfun <- function(tvec, param, Q1, Q2) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      # slice out sub-block of Q1, Q2
      ### (Verify?) Q1 is block-diagonal of size (n·iidReps) × (n·iidReps). The relevant sub-block is n×n
      Q1sub <- Q1[idx, idx, drop = FALSE]
      # Q2sub <- Q2[idx, idx, drop=FALSE]
      total <- total + cgf$K4operatorAABB(tvec[idx], param, Q1sub, Q1sub)
    }
    total
  }
  
  
  # K3K3operatorAABBCC => sum
  K3K3operatorAABBCCfun <- function(tvec, param, Q1, Q2, Q3) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      Q1sub <- Q1[idx, idx, drop = FALSE]
      # Q2sub <- Q2[idx, idx, drop=FALSE]
      # Q3sub <- Q3[idx, idx, drop=FALSE]
      total <- total + cgf$K3K3operatorAABBCC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
    }
    total
  }
  
  # K3K3operatorABCABC => sum
  K3K3operatorABCABCfun <- function(tvec, param, Q1, Q2, Q3) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    total <- 0
    for (i in seq_len(N)) {
      idx <- chunkIndices(i, bS)
      Q1sub <- Q1[idx, idx, drop = FALSE]
      # Q2sub <- Q2[idx, idx, drop=FALSE]
      # Q3sub <- Q3[idx, idx, drop=FALSE]
      total <- total + cgf$K3K3operatorABCABC(tvec[idx], param, Q1sub, Q1sub, Q1sub)
    }
    total
  }
  
  
  # ineq_constraint => concatenation
  ineq_constraintfun <- function(tvec, param) {
    n <- length(tvec)
    N <- get_nBlocks(n)
    bS <- get_blockSize(n)
    
    # call once for length
    idx_first <- chunkIndices(1, bS)
    first_val <- cgf$ineq_constraint(tvec[idx_first], param)
    L <- length(first_val)
    
    out_ <- numeric(L * N) * param[1] #### Modified to depend on `param`
    out_[1:L] <- first_val
    if (N > 1) {
      for (i in 2:N) {
        idx <- chunkIndices(i, bS)
        val <- cgf$ineq_constraint(tvec[idx], param)
        start_ <- (i - 1)*L + 1
        out_[start_:(i*L)] <- val
      }
    }
    out_
  }
  
  # for the analytic_tvec_hat_func:
  # We'll do a chunk approach if cgf$analytic_tvec_hat() is non-NULL:
  # e.g. chunk x => pass each chunk to cgf$analytic_tvec_hat => combine?
  #### Check if this doesn't make sense, (default to NULL if that's the case)
  analytic_tvec_hat_func <- NULL # If the base CGF had no valid function, just return NULL
  if (cgf$has_analytic_tvec_hat()) {
    analytic_tvec_hat_func <- function(x, param) {
      n <- length(x)
      N <- get_nBlocks(n)
      bS <- get_blockSize(n)
      out_ <- numeric(n) * param[1]
      for (i in seq_len(N)) {
        idx <- chunkIndices(i, bS)
        out_[idx] <- cgf$analytic_tvec_hat(x[idx], param)
      }
      out_
    }
  }
  
  
  # ------------------------------------------------------------------
  # Construct a character label for op_name
  # ------------------------------------------------------------------
  if (useBoth) {
    op_label <- sprintf("iidReplicatesCGF(iidReps=%d,bS=%d)", iidReps, block_size)
  } else if (onlyIID) {
    op_label <- sprintf("iidReplicatesCGF(iidReps=%d)", iidReps)
  } else { # onlyBlock
    op_label <- sprintf("iidReplicatesCGF(bS=%d)", block_size)
  }
  
  
  # ------------------------------------------------------------------
  # Build the new CGF object
  # ------------------------------------------------------------------
  createCGF(
    K = Kfun, 
    K1 = K1fun, 
    K2 = K2fun, 
    K3operator = K3operatorfun, 
    K4operator = K4operatorfun, 
    ineq_constraint = ineq_constraintfun,
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
    op_name = c(cgf$call_history, op_label)
  )

}



















#' @title Replicate a CGF object over multiple i.i.d. blocks and/or with a fixed block size
#'
#' @description
#' Extends a given `CGF` object to handle multiple i.i.d. blocks. You can specify:
#' 
#' - \code{iidReps} only,
#' - \code{block_size} only,
#' - or both \code{iidReps} and \code{block_size}.
#'
#' @param cgf A `CGF` object. 
#' @param iidReps Either \code{NULL} or a positive integer specifying the number of i.i.d. blocks.
#' @param block_size Either \code{NULL} or a positive integer specifying the size of each block.
#'
#' @return A new CGF object that operates on length-\eqn{m d} input vectors, chunked
#'   into \eqn{m} blocks, each block of length \eqn{d}.
#' 
#' @details
#' \itemize{
#'   \item If \code{iidReps} is provided (without \code{block_size}), the length of the input vector is split evenly into \code{iidReps} blocks.
#'   \item If \code{block_size} is provided (without \code{iidReps}), the input vector length determines how many blocks there are.
#'   \item If both are provided, the input vector length must be \code{block_size * iidReps}, creating exactly \code{iidReps} blocks, each of length \code{block_size}.
#' }
#'
#' @examples
#' # Suppose cgf is dimension=3, and we want 5 blocks => dimension=15
#' #   aggregator <- iidReplicatesCGF(cgf, iidReps=5)
#'
#' @export
iidReplicatesCGF <- function(cgf, iidReps = NULL, block_size = NULL) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object.")
  if (is.null(iidReps) && is.null(block_size)) stop("At least one of 'iidReps' or 'block_size' must be non-NULL.")
  if (!is.null(iidReps)) {
    if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1) {
      stop("'iidReps' must be a positive integer.")
    }
  }
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1 || block_size < 1) {
      stop("'block_size' must be a positive integer.")
    }
  }
  
  
  .iidReplicatesCGF_internal(
    cgf        = cgf, 
    iidReps    = iidReps,
    block_size = block_size
  )

}






