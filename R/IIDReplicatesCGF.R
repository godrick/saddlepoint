# R/IIDReplicatesCGF.R
# Objects: iidReplicatesCGF
# Areas to review are marked with `####`



#------------------------------------------------------------
# Helper functions 
#------------------------------------------------------------
# a helper for chunking
chunkIndices <- function(i, block_size) {
  seq.int((i - 1)*block_size + 1, i*block_size)
}


#' @title Replicate a CGF over Multiple IID Blocks
#'
#' @description
#' This function takes any valid CGF object and extends it to handle
#' multiple IID blocks. Given \code{iidReps}, the input \code{tvec} is treated
#' as \code{iidReps} separate blocks, each of size \eqn{block_size = length(tvec)/\mathrm{iidReps}}.
#'
#' For each operator (\eqn{K, K1, K2}, etc.), \code{iidReplicatesCGF()}:
#' \itemize{
#'   \item Splits \code{tvec} (and any corresponding vectors \code{x,y,v1,...}) into sub-blocks.
#'   \item Calls the original CGF operator on each sub-block.
#'   \item Combines results by summation (for \eqn{K}, \eqn{K3}, \eqn{K4}), or by
#'         forming block-diagonal (for \eqn{K2} and similar).
#' }
#'
#' Effectively, if the original CGF is compatible with dimension \eqn{d}, and \code{iidReps=m},
#' the resulting CGF will be compatible with dimension \eqn{m d}, chunked into \eqn{m} blocks each
#' of size \eqn{d}.
#'
#' @param cgf A CGF object. Must be dimension \eqn{d}.
#' @param iidReps A positive integer \eqn{m} specifying how many IID blocks to create.
#'   Must be at least 2 (since \code{iidReps=1} would just be the original CGF).
#' @param ... Additional arguments passed to \code{\link{createCGF}} in case
#'   you want to override any methods.
#'
#' @return A new CGF object that operates on length-\eqn{m d} input vectors, chunked
#'   into \eqn{m} blocks, each block of length \eqn{d}.
#'
#' @details
#' Typically, you do not call `iidReplicatesCGF()` directly; instead, you set
#' \code{iidReps > 1} in a user-facing function like
#' \code{\link{linearlyMappedCGF}}, which internally calls `iidReplicatesCGF()`.
#'
#' @examples
#' # Suppose cgf is dimension=3, and we want 5 blocks => dimension=15
#' #   aggregator <- iidReplicatesCGF(cgf, iidReps=5)
#'
#' @export
iidReplicatesCGF <- function(cgf, iidReps, ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object.")
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 2) stop("'iidReps' must be an integer >= 2.")
  
  # fetch some private methods from the base CGF
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  neg_ll <- cgf$.get_private_method("neg_ll")
  func_T <- cgf$.get_private_method("func_T")
  
  
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
  if ( cgf$has_analytic_tvec_hat() ) {
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






