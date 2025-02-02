# R/ConcatenationCGF.R

#' @title Concatenation of CGF objects.
#'
#' @description
#' Constructs a CGF object for a concatenated vector of independent random variables, 
#' with each represented by its respective CGF object and a corresponding dimension. 
#' This means CGF \code{cgf_list[[i]]} sees \code{block_sizes[i]} consecutive components of \code{tvec}.
#' It facilitates the modeling of a concatenated random vector, 
#' where each segment of the vector is specified by its CGF object and a corresponding dimension.
#'
#' @param cgf_list A list of CGF objects, each inheriting from class \code{"CGF"}.
#' @param block_sizes An integer vector, same length as \code{cgf_list}, specifying how many
#'   components of \code{tvec} each CGF expects.
#' @param iidReps A positive integer. If greater than 1, replicates the resulting CGF \code{iidReps} times. Defaults to 1.
#' @param adaptor A function(\code{theta}) -> numeric vector, transforming model parameter vector
#'   \code{theta} into what the resulting CGF object expects. By default, the identity function.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, 
#'   or to \code{\link{iidReplicatesCGF}} if \code{iidReps>1}.
#'
#' 
#' @seealso \code{\link{iidReplicatesCGF}}, \code{\link{sumOfIndependentCGF}}
#'
#' @return A CGF object.
#'
#' @examples
#' \dontrun{
#' # Suppose we have two CGFs with dimension 3 and 2, respectively:
#' # CGF1 sees tvec of length=3, CGF2 sees tvec of length=2
#' # We want to form a blockwise CGF with dimension=3+2=5.
#'
#' cgf_obj1 <- someCGF_dim3
#' cgf_obj2 <- someCGF_dim2
#' # Then:
#' concat_cgf <- concatenateBlockwiseCGF(
#'   cgf_list = list(cgf_obj1, cgf_obj2),
#'   block_sizes = c(3,2)
#' )
#'
#' # Evaluate K at tvec of length=5
#' tvec_example <- c(0.1, -0.2, 0.05, 0.5, 0.4)
#' concat_cgf$K(tvec_example, param=c(...))
#'
#' # If we want i.i.d. replicates of that block structure => dimension= 5*iidReps
#' concat_cgf2 <- concatenateBlockwiseCGF(
#'   cgf_list = list(cgf_obj1, cgf_obj2),
#'   block_sizes = c(3,2),
#'   iidReps = 2
#' )
#' # Now tvec must be length=10, etc.
#' }
#' @export
concatenationCGF <- function(cgf_list, 
                             block_sizes,
                             iidReps = 1,
                             adaptor = function(x) x,
                             ...) {
  # checks
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if (length(block_sizes) != length(cgf_list)) stop("'block_sizes' must have the same length as 'cgf_list'.")
  if (any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE))) stop("Every element of 'cgf_list' must be an object of class 'CGF'.")
  if (!all(block_sizes == as.integer(block_sizes)) || any(block_sizes <= 0)) stop("'block_sizes' must be a vector of positive integers.")
  if (!is.numeric(iidReps) || length(iidReps)!=1 || iidReps<1) stop("'iidReps' must be a positive integer.")
  
  total_dim <- sum(block_sizes)
  
  # K => sum of child K, but each sees sub-block
  Kfun <- function(tvec, param) {
    # param is global, each CGF in cgf_list knows what to do with it
    # tvec must have length= total_dim
    if (length(tvec)!= total_dim) stop(sprintf("`tvec` has length %d, expected %d from block_sizes sum.", length(tvec), total_dim))
    
    total <- 0
    current_start <- 1
    for (i in seq_along(cgf_list)) {
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
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
      len_i <- block_sizes[i]
      idx   <- current_start:(current_start + len_i - 1)
      piece <- cgf_list[[i]]$ineq_constraint(tvec[idx], param)
      
      # Append child's constraints to out_
      if (length(piece) > 0) out_ <- c(out_, piece)
      current_start <- current_start + len_i
    }
    
    out_
  }
  
  
  # paste(vapply(cgf_list, function(cg) paste(cg$call_history, collapse="->"), ""), collapse=","),
  
  # call_history => for debugging
  combined_history <- paste0("{", 
                             paste(vapply(cgf_list, 
                                          function(x) paste0(x$call_history, collapse=" -> "), 
                                          FUN.VALUE=""), 
                                   collapse=", "), 
                             "}")
  op_name_vec <- c(combined_history, "concatenationCGF")
  param_adaptor_ <- validate_function_or_adaptor(adaptor)
  
  # create the CGF
  final_cgf <- createCGF(
    K = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    
    func_T = funcTfun,
    
    ineq_constraint = ineqfun,
    param_adaptor   = param_adaptor_,
    
    op_name = op_name_vec,
    ...
  )
  
  if (iidReps == 1) return(final_cgf)
  
  
  iidReplicatesCGF(cgf = final_cgf, iidReps = iidReps, ...)
  
}
