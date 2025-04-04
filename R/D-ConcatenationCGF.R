# R/ConcatenationCGF.R



.concatenationCGF_internal <- function(cgf_list, 
                                       component_dims,
                                       ...){
  
  total_dim <- sum(component_dims)
  
  # K => sum of child K, but each sees sub-block
  Kfun <- function(tvec, param) {
    # param is global, each CGF in cgf_list knows what to do with it
    # tvec must have length = total_dim
    if (length(tvec)!= total_dim) stop(sprintf("`tvec` has length %d, expected %d from component_dims sum.", length(tvec), total_dim))
    total <- 0
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      t_sub <- tvec[current_start:(current_start+len_i-1)]
      total <- total + cgf_list[[i]]$K(t_sub, param)
      current_start <- current_start + len_i
    }
    total
  }
  
  
  # K1 => piecewise concatenation
  K1fun <- function(tvec, param) {
    if (length(tvec)!= total_dim) stop(sprintf("`tvec` has length %d, expected %d.", length(tvec), total_dim))
    
    out <- numeric(total_dim) * param[1]  #### Modified to depend on `param` // temporary fix for RTMB error
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      t_sub <- tvec[current_start:(current_start+len_i-1)]
      out_sub <- cgf_list[[i]]$K1(t_sub, param)
      out[current_start:(current_start+len_i-1)] <- out_sub
      current_start <- current_start + len_i
    }
    out
  }
  
  # K2 => block diagonal
  K2fun <- function(tvec, param) {
    if (length(tvec)!= total_dim) stop(sprintf("`tvec` length mismatch: got %d, expected %d", length(tvec), total_dim))
    
    accum <- matrix(0, nrow=total_dim, ncol=total_dim) * param[1]  #### Modified to depend on `param` // temporary fix for RTMB error
    ##### possibly a sparse matrix??; But Matrix::determinant() fails with adsparse matrices.
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      t_sub <- tvec[current_start:(current_start+len_i-1)]
      k2_sub <- cgf_list[[i]]$K2(t_sub, param)
      # place it in accum's block
      accum[current_start:(current_start+len_i-1), current_start:(current_start+len_i-1)] <- k2_sub
      current_start <- current_start + len_i
    }
    accum
  }
  
  # K3operator => sum of child K3operator
  K3opfun <- function(tvec, param, v1, v2, v3) {
    if (length(tvec)!= total_dim || length(v1)!= total_dim ||
        length(v2)!= total_dim || length(v3)!= total_dim) {
      stop("dimension mismatch in K3operator arguments.")
    }
    total <- 0
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      idx <- current_start:(current_start+len_i-1)
      total <- total + cgf_list[[i]]$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
      current_start <- current_start + len_i
    }
    total
  }
  
  # K4operator => sum across sub-block operators
  K4opfun <- function(tvec, param, v1, v2, v3, v4) {
    if (length(tvec)!= total_dim || length(v1)!= total_dim ||
        length(v2)!= total_dim || length(v3)!= total_dim || length(v4)!= total_dim) {
      stop("dimension mismatch in K4operator arguments.")
    }
    total <- 0
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      idx <- current_start:(current_start+len_i-1)
      total <- total + cgf_list[[i]]$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
      current_start <- current_start + len_i
    }
    total
  }
  
  
  # func_T => sum of child func_T
  func_T_list <- lapply(cgf_list, function(cg) cg$.get_private_method("func_T"))
  funcTfun <- function(tvec, param) {
    if (length(tvec)!= total_dim) {
      stop(sprintf("`tvec` length mismatch in func_T: got %d, expected %d", 
                   length(tvec), total_dim))
    }
    total_res <- 0
    current_start <- 1
    for (i in seq_along(func_T_list)) {
      len_i <- component_dims[i]
      idx <- current_start:(current_start + len_i - 1)
      # Typically, child$func_T(t_sub, param) is a scalar
      total_res <- total_res + func_T_list[[i]](tvec[idx], param)
      current_start <- current_start + len_i
    }
    total_res
  }
  
  # ineq_constraint => concatenation
  ineqfun <- function(tvec, param) {
    if (length(tvec) != total_dim) stop(sprintf("`tvec` length mismatch in ineq_constraint: got %d, expected %d", length(tvec), total_dim))
    
    # We'll build final constraints in a single pass, 
    # appending for each child (no repeated calls).
    ##### This is risky, as the length of the output is unknown; please test/check 
    out_ <- numeric(0)*param[1]  #### Modified to depend on `param` // temporary fix for RTMB error
    
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- component_dims[i]
      idx   <- current_start:(current_start + len_i - 1)
      piece <- cgf_list[[i]]$ineq_constraint(tvec[idx], param)
      
      # Append child's constraints to out_
      if (length(piece) > 0) out_ <- c(out_, piece)
      current_start <- current_start + len_i
    }
    
    out_
  }
  
  # call_history => for debugging
  combined_history <- paste0("{", 
                             paste(vapply(cgf_list, 
                                          function(x) paste0(x$call_history, collapse=" -> "), 
                                          FUN.VALUE=""), 
                                   collapse=", "), 
                             "}")
  op_name_vec <- c(combined_history, "concatenationCGF")
  
  createCGF(
    K = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    
    func_T = funcTfun,
    ineq_constraint = ineqfun,
    op_name = op_name_vec,
    ...
  )
  
  
  
}



#' @title Concatenation of CGF objects.
#'
#' @description
#' Constructs a CGF object for a concatenated vector of independent random variables, 
#' with each represented by its respective CGF object and a corresponding dimension. 
#' This means CGF \code{cgf_list[[i]]} sees \code{component_dims[i]} consecutive components of \code{tvec}.
#' It facilitates the modeling of a concatenated random vector, 
#' where each segment of the vector is specified by its CGF object and a corresponding dimension.
#' 
#'
#' @param cgf_list A non-empty list of `CGF` objects.
#' @param component_dims An integer vector (same length as \code{cgf_list}) specifying the length of \code{tvec} compatible with each `CGF`. Defaults to 1 for each `CGF`.
#' @param block_size Optional integer specifying the block size for i.i.d. replication of the resulting `CGF`. Defaults to \code{NULL}.
#' @param iidReps Optional integer specifying the number of i.i.d. blocks that will be expected the resulting `CGF`. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to the `CGF` creation function.
#' 
#' @details
#' - The dimension of \code{tvec} is \code{sum(component_dims)}.
#' - If \code{iidReps} or \code{block_size} is specified, the result is automatically wrapped
#'   via \code{iidReplicatesCGF()}.
#'   
#' 
#' @seealso \code{\link{iidReplicatesCGF}}, \code{\link{sumOfIndependentCGF}}
#'
#' @return A `CGF` object.
#'
#' @examples
#' \dontrun{
#' # 
#' }
#' @export
concatenationCGF <- function(cgf_list, 
                             component_dims = 1,
                             block_size = NULL, 
                             iidReps = NULL,
                             ...) {
  # checks
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if (any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE))) stop("Every element of 'cgf_list' must be an object of class 'CGF'.")
  
  if (!all(component_dims == as.integer(component_dims)) || any(component_dims <= 0)) {
    stop("'component_dims' must be positive integers.")
  } else if (length(component_dims) == 1) {
    component_dims <- rep(component_dims, length(cgf_list))
  } else if (length(component_dims) != length(cgf_list)) {
    stop("'component_dims' must have the same length as 'cgf_list'.")
  }
  
  base_cgf <- .concatenationCGF_internal(cgf_list = cgf_list, component_dims = component_dims, ...)

  if (is.null(block_size) && is.null(iidReps)) return(base_cgf)
  if (!is.null(iidReps) && iidReps == 1) return(base_cgf)
  
  iidReplicatesCGF(cgf = base_cgf, iidReps = iidReps, block_size = block_size)

}
