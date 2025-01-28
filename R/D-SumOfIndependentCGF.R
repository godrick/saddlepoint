# ------------------------------------------------------------
#  R/SumOfIndependentCGF.R
#  Function: sumOfIndependentCGF, .sumOfIndependentCGF_internal
#  Purpose: Create a CGF object for the the sum of independent rvs,
#           with optional replication for i.i.d. observations.
# ------------------------------------------------------------




.sumOfIndependentCGF_internal <- function(cgf_list, ...){
  #-------------------------------------
  # First, if both iidReps and block_size are NULL => no replication
  #-------------------------------------
  
  #  Combine each method: K, K1, K2, etc.
  #  We do a simple for loop accumulation.
  #  For scalar results (like K, K3operator, etc.), we sum them.
  #  For vector results (K1, etc.), we do elementwise sum.
  #  For matrix results (K2, etc.), we do matrix sum.
  #  For constraints, we concatenate them.
  
  ## K
  Kfun <- function(tvec, param) {
    total <- 0
    for (cg in cgf_list) {
      total <- total + cg$K(tvec, param)
    }
    total
  }
  
  ## K1
  K1fun <- function(tvec, param) {
    out <- numeric(length(tvec))
    for (cg in cgf_list) {
      out <- out + cg$K1(tvec, param)
    }
    out
  }
  
  ## K2
  K2fun <- function(tvec, param) {
    dim_ <- length(tvec)
    accum <- matrix(0, nrow = dim_, ncol = dim_) ##### possibly a sparse matrix??; But Matrix::determinant() fails with adsparse matrices.
    for (cg in cgf_list) {
      accum <- accum + cg$K2(tvec, param)
    }
    accum
  }
  
  ## K3operator
  K3opfun <- function(tvec, param, v1, v2, v3) {
    total <- 0
    for (cg in cgf_list) {
      total <- total + cg$K3operator(tvec, param, v1, v2, v3)
    }
    total
  }
  
  ## K4operator
  K4opfun <- function(tvec, param, v1, v2, v3, v4) {
    total <- 0
    for (cg in cgf_list) {
      total <- total + cg$K4operator(tvec, param, v1, v2, v3, v4)
    }
    total
  }
  
  
  tilitingEx_list <- lapply(cgf_list, function(x) x$.get_private_method("tilting_exponent"))
  # negll_list <- lapply(cgf_list, function(x) x$.get_private_method("neg_ll"))
  # func_T_list <- lapply(cgf_list, function(x) x$.get_private_method("func_T"))
  
  tiltingfun <- function(tvec, param) {
    total <- 0
    for (cg in tilitingEx_list) {
      total <- total + cg(tvec, param)
    }
    total
  }
  
  
  
  
  ## K2operator
  K2opfun <- function(tvec, param, x, y) {
    total <- 0
    for (cg in cgf_list) {
      total <- total + cg$K2operator(tvec, param, x, y)
    }
    total
  }
  
  ## K2operatorAK2AT
  K2opAK2ATfun <- function(tvec, param, B) {
    # We'll sum up the resulting matrices
    dim_ <- nrow(B)
    accum <- Matrix(0, nrow = dim_, ncol = dim_) # possibly a sparse matrix??
    for (cg in cgf_list) {
      accum <- accum + cg$K2operatorAK2AT(tvec, param, B)
    }
    accum
  }
  
  ## K4operatorAABB
  K4AABBfun <- function(tvec, param, Q1, Q2) {
    total <- 0
    for (cg in cgf_list) {
      total <- total + cg$K4operatorAABB(tvec, param, Q1, Q2)
    }
    total
  }
  
  
  
  
  ## ineq_constraint
  # We'll just concatenate them
  ineqfun <- function(tvec, param) {
    # Calculate the total size needed for the result
    # total_size <- 0
    # for (cg in cgf_list) {
    #   total_size <- total_size + length(cg$ineq_constraint(tvec, param))
    # }
    total_size <- sum(sapply(cgf_list, function(cg) length(cg$ineq_constraint(tvec, param))))
    out_ <- numeric(total_size)*param[1] #### A temporary fix to avoid RTMB error. See IIDReplicatesCGF.R for more details.
    if(total_size > 0){
      start_idx <- 1
      for (cg in cgf_list) {
        piece <- cg$ineq_constraint(tvec, param)
        if(length(piece) > 0){
          end_idx <- start_idx + length(piece) - 1
          out_[start_idx:end_idx] <- piece
          start_idx <- end_idx + 1
        }
      }
    }
    
    out_
  }
  
  
  # Build a combined 'call_history' vector:
  # We'll concatenate the call_history from each CGF in the final list
  combined_history <- paste0("[", paste(sapply(cgf_list, function(cg) cg$call_history), collapse = ", "), "]" )
  op_name_vec <- c(combined_history, "sumOfIndependentCGF")
  
  createCGF(
    K = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    
    
    tilting_exponent = tiltingfun,
    ineq_constraint = ineqfun,
    
    K2operator = K2opfun,
    K2operatorAK2AT = K2opAK2ATfun,
    
    # func_T = func_Tfun,
    K4operatorAABB = K4AABBfun,
    # K3K3operatorAABBCC = K3K3AABBCCfun,
    # K3K3operatorABCABC = K3K3ABCABCfun,
    
    
    op_name = op_name_vec,
    ... # pass in any further overrides
  )
}





#' @title CGF Object for the sum of independent random variables
#'
#' @description
#' Constructs a new CGF object representing the sum of multiple independent random variables.
#'
#'
#' @param cgf_list A non-empty list of CGF objects. Each of class \code{"CGF"}.
#' @param block_size Either \code{NULL} or a positive integer specifying the block size for replication.
#'   Default is \code{NULL}.
#' @param iidReps Either \code{NULL} or a positive integer specifying how many i.i.d. blocks 
#'   to expect. Default is \code{NULL}.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, the `CGF` object creation function.
#' 
#' @seealso \code{\link{adaptCGF}} for an example which also adapts the CGF object.
#' 
#' @return A new `CGF` object.
#' @export
sumOfIndependentCGF <- function(cgf_list, block_size = NULL, iidReps = NULL, ...) {
  # Basic validations
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if ( any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE)) ) stop("Every element of 'cgf_list' must be of class 'CGF'." )
  
  cgf_res <- .sumOfIndependentCGF_internal(cgf_list, ...)
  
  if (is.null(block_size) && is.null(iidReps)) return(cgf_res)
  if (!is.null(iidReps) && iidReps == 1) return(cgf_res)
  iidReplicatesCGF(cgf = cgf_res, iidReps = iidReps, block_size = block_size)
}


