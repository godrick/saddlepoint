# ------------------------------------------------------------------
#  R/D-RandomlyStoppedSumCGF.R
#  Main function: randomlyStoppedSumCGF (exported)
#  Internal : .randomlyStoppedSumCGF_internal
#  Purpose: Create a CGF object for Y = sum_{i=1}^N X_i,
#           where N ~ count_cgf, and X_i ~ summand_cgf (i.i.d.),
#           with N independent of the X_i's.
# ------------------------------------------------------------------



.randomlyStoppedSumCGF_internal <- function(count_cgf, summand_cgf, ...){
  # ------------------------------------------------------------------
  # Build the methods for the new CGF
  # ------------------------------------------------------------------
  
  # Utility: turn the scalar K_X(t) into a single-element "vector" for count_cgf
  # so that count_cgf$K() sees tvec dimension=1
  # We also consistently treat K1_N(tCount) as a 1-element vector or 1x1 matrix
  # that we interpret as a scalar.
  
  # K_Y = K_N( K_X(tvec, params), params)
  Kfun <- function(tvec, param) {
    # summand_cgf$K(tvec, param) is a scalar => call it scalar_Kx
    scalar_Kx <- summand_cgf$K(tvec, param)
    count_cgf$K(scalar_Kx, param)
  }
  
  # K1_Y = K1_N( K_X(tvec, params), params) * K1_X(tvec, params)
  K1fun <- function(tvec, param) {
    scalar_Kx <- summand_cgf$K(tvec, param)
    
    K1N_vec <- count_cgf$K1(scalar_Kx, param)  # length=1 numeric
    K1X_vec <- summand_cgf$K1(tvec, param)     # same dimension as tvec
    # multiply K1N by each entry in K1X
    K1N_vec * K1X_vec
  }
  
  # K2_Y = (K1_X) (**) (K1_X)^T + (*) K2_X
  # *  => K1_N(K_X(tvec, params), params) // this will be a scalar
  # ** => K2_N(K_X(tvec, params), params) // this is also a scalar
  K2fun <- function(tvec, param) {
    scalar_Kx <- summand_cgf$K(tvec, param)
    K1N_vec <- count_cgf$K1(scalar_Kx, param)  # length=1
    K2N_mat <- count_cgf$K2(scalar_Kx, param)  # 1x1
    K1X_vec <- summand_cgf$K1(tvec, param)     # length(tvec)
    K2X_mat <- summand_cgf$K2(tvec, param)     # dimension= length(tvec) x length(tvec)
    
    (K1X_vec %*% K2N_mat %*% t(K1X_vec)) + (K1N_vec * K2X_mat)
  }
  
  K3opfun <- function(tvec, param, w1, w2, w3) {
    scalar_Kx <- summand_cgf$K(tvec, param)
    K1N_vec <- count_cgf$K1(scalar_Kx, param) # length=1
    K1X_vec <- summand_cgf$K1(tvec, param)    # length(tvec)
    
    K1N_vec * summand_cgf$K3operator(tvec, param, w1, w2, w3) +
      count_cgf$K2operator(scalar_Kx, param, sum(w1*K1X_vec), summand_cgf$K2operator(tvec, param, w2, w3) )  +
      count_cgf$K2operator(scalar_Kx, param, sum(w2*K1X_vec), summand_cgf$K2operator(tvec, param, w1, w3) )  +
      count_cgf$K2operator(scalar_Kx, param, sum(w3*K1X_vec), summand_cgf$K2operator(tvec, param, w1, w2) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w1*K1X_vec), sum(w2*K1X_vec), sum(w3*K1X_vec) )
  }
  
  
  
  K4opfun <- function(tvec, param, w1, w2, w3, w4) {
    scalar_Kx <- summand_cgf$K(tvec, param)
    K1N_vec <- count_cgf$K1(scalar_Kx, param)  # length=1
    K1X_vec <- summand_cgf$K1(tvec, param)     # length(tvec)
    
    K1N_vec*summand_cgf$K4operator(tvec, param, w1, w2, w3, w4) +
      count_cgf$K2operator(scalar_Kx, param, sum(w1*K1X_vec), summand_cgf$K3operator(tvec, param, w2, w3, w4) )  +
      count_cgf$K2operator(scalar_Kx, param, sum(w2*K1X_vec), summand_cgf$K3operator(tvec, param, w1, w3, w4) )  +
      count_cgf$K2operator(scalar_Kx, param, sum(w3*K1X_vec), summand_cgf$K3operator(tvec, param, w1, w2, w4) )  +
      count_cgf$K2operator(scalar_Kx, param, sum(w4*K1X_vec), summand_cgf$K3operator(tvec, param, w1, w2, w3) )  +
      count_cgf$K2operator(scalar_Kx, param, summand_cgf$K2operator(tvec, param, w1, w2), summand_cgf$K2operator(tvec, param, w3, w4) )  +
      count_cgf$K2operator(scalar_Kx, param, summand_cgf$K2operator(tvec, param, w1, w4), summand_cgf$K2operator(tvec, param, w2, w3) )  +
      count_cgf$K2operator(scalar_Kx, param, summand_cgf$K2operator(tvec, param, w1, w3), summand_cgf$K2operator(tvec, param, w2, w4) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w1*K1X_vec), sum(w2*K1X_vec), summand_cgf$K2operator(tvec, param, w3, w4) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w1*K1X_vec), sum(w3*K1X_vec), summand_cgf$K2operator(tvec, param, w2, w4) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w1*K1X_vec), sum(w4*K1X_vec), summand_cgf$K2operator(tvec, param, w2, w3) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w2*K1X_vec), sum(w3*K1X_vec), summand_cgf$K2operator(tvec, param, w1, w4) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w2*K1X_vec), sum(w4*K1X_vec), summand_cgf$K2operator(tvec, param, w1, w3) )  +
      count_cgf$K3operator(scalar_Kx, param, sum(w3*K1X_vec), sum(w4*K1X_vec), summand_cgf$K2operator(tvec, param, w1, w2) )  +
      count_cgf$K4operator(scalar_Kx, param, sum(w1*K1X_vec), sum(w2*K1X_vec), sum(w3*K1X_vec), sum(w4*K1X_vec) )
    
  }
  
  
  K2opfun <- function(tvec, param, x, y) {
    scalar_Kx <- summand_cgf$K(tvec, param)
    K1N_vec <- count_cgf$K1(scalar_Kx, param)
    K1X_vec <- summand_cgf$K1(tvec, param)
    
    K1N_vec * summand_cgf$K2operator(tvec, param, t(x), y) +
      count_cgf$K2operator(scalar_Kx, param, sum(K1X_vec*x), sum(K1X_vec*y) )
  }
  
  
  ineqfun <- function(tvec, param) {
    # summand_ineq is dimension= length(tvec)
    summand_ineq <- summand_cgf$ineq_constraint(tvec, param)
    # count_ineq is dimension for 1D?
    scalar_Kx <- summand_cgf$K(tvec, param)
    count_ineq <- count_cgf$ineq_constraint(scalar_Kx, param)
    
    c(count_ineq, summand_ineq)
  }
  
  # ------------------------------------------------------------------
  # Create the new CGF object
  # ------------------------------------------------------------------
  # A small call_history label:
  combined_history <- paste0(
    "count: ", count_cgf$call_history, "\n",
    "summand: ", summand_cgf$call_history
  )
  op_name_vec <- c(combined_history, "randomlyStoppedSumCGF")
  
  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    K2operator = K2opfun,
    ineq_constraint = ineqfun,
    op_name = op_name_vec,
    ...
  )
  
}




#' @title CGF Object for a Randomly-Stopped Sum
#'
#' @description 
#' Constructs the CGF for a randomly stopped sum 
#' \eqn{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}{Y = sum(X[1], X[2], ..., X[N])},
#' where \eqn{N} is a non-negative integer-valued random variable with CGF object 
#' specified by \code{count_cgf}, and each \eqn{\tilde{X}_i} is an i.i.d. random 
#' vector with CGF specified by \code{summand_cgf}. 
#' \eqn{N} is independent of \eqn{\tilde{X}_i}.
#' 
#' Optionally, the resulting CGF object can be replicated in i.i.d. blocks if 
#' \code{iidReps} or \code{block_size} (or both) are specified.
#' 
#' 
#' @param count_cgf A `CGF` object corresponding to the distribution of \eqn{N} 
#'   in the sum \eqn{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}.
#' @param summand_cgf A `CGF` object corresponding to the summand random variables 
#'   \eqn{\tilde{X}_i} in the sum.
#' @param block_size (Optional) A positive integer specifying the size of each 
#'   block i.i.d. replication; or \code{NULL} (the default) to skip block-size replication.
#' @param iidReps (Optional) A positive integer specifying how many i.i.d. blocks 
#'   to expect; or \code{NULL} (the default) to skip i.i.d. replication.
#' @param ... Additional named arguments passed to CGF creation functions.
#'
#' @details
#' **Setup**: 
#' \deqn{{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}}
#' The random variable \eqn{\tilde{Y}} is computed as the sum of \eqn{N} i.i.d copies of the random vector \eqn{\tilde{X}_i}. Each \eqn{\tilde{X}_i} is a random vector of arbitrary dimension \eqn{d}.
#' For example, when \eqn{d=3}:
#' \deqn{ \tilde{X_1} = (X_{11}, X_{12}, X_{13}) }
#' \deqn{ \tilde{X_2} = (X_{21}, X_{22}, X_{23}) }
#' \deqn{ ... }
#' \deqn{ \tilde{X_N} = (X_{N1}, X_{N2}, X_{N3}) }
#' Here, \eqn{X_{i1}, X_{i2}, X_{i3}} are scalars but the vectors \eqn{\tilde{X_i}} are i.i.d copies.
#'
#' **Optional Replication**:
#' If \code{block_size} or \code{iidReps} (or both) are specified, the function 
#' automatically calls \code{\link{iidReplicatesCGF}} on the resulting CGF.  
#' You may specify:
#' \enumerate{
#'   \item \code{iidReps} only,
#'   \item \code{block_size} only,
#'   \item **both** \code{iidReps} and \code{block_size}.
#' }
#' If both are provided, the input dimension of the final CGF must be 
#' \code{block_size * iidReps}.
#'
#'
#' @return A `CGF` object for the distribution of \eqn{\sum_{i=1}^N \tilde{X}_i}.
#' @export
#'
#' @examples
#' \dontrun{
#' # For example, let count_cgf be the CGF of N ~ Poisson(lambdaN),
#' # and summand_cgf be the CGF of X ~ Gamma(shape, rate).
#' # Then Y = sum_{i=1..N} X_i is a compound Poisson-Gamma distribution.
#'
#' cgf_Y <- randomlyStoppedSumCGF(count_cgf = PoissonCGF, summand_cgf = GammaCGF)
#'
#'
#' }
randomlyStoppedSumCGF <- function(count_cgf,
                                  summand_cgf,
                                  block_size = NULL,
                                  iidReps    = NULL,
                                  ...) 
{
  if (!inherits(count_cgf, "CGF")) stop("'count_cgf' must be a CGF object.")
  if (!inherits(summand_cgf, "CGF")) stop("'summand_cgf' must be a CGF object.")
  out_cgf <- .randomlyStoppedSumCGF_internal(count_cgf, summand_cgf, ...)
  
  if (is.null(block_size) && is.null(iidReps)) return(out_cgf)
  if (!is.null(iidReps) && iidReps == 1) return(out_cgf)
  iidReplicatesCGF(cgf = out_cgf, iidReps = iidReps, block_size = block_size)
}
