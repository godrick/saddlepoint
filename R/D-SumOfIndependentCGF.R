# ------------------------------------------------------------
#  R/SumOfIndependentCGF.R
#  Function: sumOfIndependentCGF
#  Purpose: Create a CGF object for the the sum of independent rvs,
#           with optional replication for i.i.d. observations.
# ------------------------------------------------------------



#' @title CGF Object for the sum of independent random variables
#'
#' @description
#' Constructs a new CGF object representing the sum of multiple independent random variables.
#' An optional adaptor function can be provided to transform the global parameter vector into the specific parameters required by the resulting CGF.
#'
#'
#' @param cgf_list A non-empty list of CGF objects. Each of class \code{"CGF"}.
#' @param block_size Either \code{NULL} or a positive integer specifying the block size for replication.
#'   Default is \code{NULL}.
#' @param iidReps Either \code{NULL} or a positive integer specifying how many i.i.d. blocks 
#'   to expect. Default is \code{NULL}.
#' @param adaptor A function to transform the global parameter vector \code{theta} into the parameter vector expected by the summed CGF.
#'   It should accept a numeric vector \code{theta} and return a transformed numeric vector.
#'   The default behavior is an identity function.
#' @param ... Additional named arguments passed to \code{\link{createCGF}} or \code{\link{iidReplicatesCGF}}.
#'
#'
#' @details
#' - **Adaptor Function**: Transforms the global parameter vector \code{theta} to match the structure required by the resulting CGF.
#'   This is essential when the summed CGF expects parameters in a different order or combination than the original model parameters.
#'   For example, if summing two r.vs where one depends on \code{lambda1 = alpha + beta} and the other on \code{lambda2 = beta + gamma},
#'   the adaptor function can will help in mapping a three-dimensional model parameter \code{theta = c(alpha, beta, gamma)} to \code{c(lambda1, lambda2)}.
#'   See the example below.
#'   
#' - **Replication for i.i.d. Observations**: This is handles in one of two ways
#' \enumerate{
#'   \item By `iidReps`: if \code{iidReps} is a positive integer.
#'   \item By `block_size`: if \code{block_size} is a positive integer.
#' }
#' If both \code{iidReps} and \code{block_size} are \code{NULL}, no replication is performed;
#'
#'
#' @examples
#' \dontrun{
#' # Example: Summing two non-identical Poisson r.vs (we use an adaptor function)
#'
#' # Scenario:
#' # Y1 ~ Poisson(alpha)
#' # Y2 ~ Poisson(2*alpha)
#' # Y = Y1 + Y2 ~ Poisson(3*alpha)
#' # Model Parameters: alpha
#' # CGF for Y expects: distribution_params = c(lambda1, lambda2) = c(alpha, 2*alpaha)
#' 
#'
#' # First, define an adaptor function to map alpha
#' # to distribution_params = c(alpha, 2*alpha) for the CGF of Y
#' mapThetaToDistParams <- function(theta) c(theta, 2*theta)
#' 
#'
#' # Next, create individual Poisson CGFs with separate lambda parameters
#' # PoissonModelCGF expects a single lambda parameter,
#' # we use its adaptor argument to ensure that each CGF object receives the correct parameter from distribution_params
#' K_Y1 <- PoissonModelCGF(lambda = function(x) x)       # For Y1: lambda1 = alpha
#' K_Y2 <- PoissonModelCGF(lambda = function(x) 2*x )    # For Y2: lambda2 = 2*alpha
#'
#' # Then the CGF for the summed rvs using sumOfIndependentCGF with the adaptor
#' sum_cgf <- sumOfIndependentCGF(
#'   cgf_list = list(K_Y1, K_Y2),
#'   adaptor = mapThetaToDistParams,
#'   iidReps = 1  # No replication in this example
#' )
#'
#' # # (Optional) Replicate the CGF for multiple i.i.d. observations
#' # # For example, replicate the resulting CGF to be compatible with 5 iid copies of Y.
#' # sum_cgf <- sumOfIndependentCGF(
#' #   cgf_list = list(K_Y1, K_Y2),
#' #   adaptor  = mapThetaToDistParams,
#' #   iidReps  = 5  
#' # )
#' # # General testing 
#' # # a sample of counts from two Poisson-distributed variables
#' # Y <- c(24, 25, 13, 30, 33)
#' # 
#' # # Using `sum_cgf`
#' # mle_sum <- find.saddlepoint.MLE(
#' #   observed.data = Y,
#' #   cgf = sum_cgf,
#' #   starting.theta = 1
#' # )$MLEs.theta
#' # 
#' # # Using a single Poisson CGF object (`pois_cgf`)
#' # # This function maps the parameter `alpha` to the Poisson rate `lambda = 3 * alpha`
#' # pois_adaptor <- function(alpha) 3 * alpha
#' # pois_cgf <- PoissonModelCGF(lambda = pois_adaptor, iidReps = 5 )
#' # mle_pois <- find.saddlepoint.MLE(
#' #   observed.data = Y,
#' #   cgf = pois_cgf,
#' #   starting.theta = 1
#' # )$MLEs.theta
#' # # Both MLEs (`mle_sum` and `mle_pois`) should yield the same or very similar results
#' }
#'
#' @return A CGF object
#' @export
sumOfIndependentCGF <- function(cgf_list, block_size = NULL, iidReps = NULL, adaptor = NULL, ...) {
  # Basic validations
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if ( any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE)) ) stop("Every element of 'cgf_list' must be of class 'CGF'." )
  
  # Only one of (iidReps, block_size) can be set:
  if (!is.null(iidReps) && !is.null(block_size)) stop("Please specify only one of 'iidReps' or 'block_size', not both.")
  
  # if iidReps is set:
  if (!is.null(iidReps)) {
    if (!is.numeric(iidReps) || length(iidReps) != 1 ||
        iidReps < 1 || iidReps != as.integer(iidReps)) {
      stop("'iidReps' must be NULL or a positive integer.")
    }
  }
  
  # if block_size is set:
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1 ||
        block_size < 1 || block_size != as.integer(block_size)) {
      stop("'block_size' must be NULL or a positive integer.")
    }
  }
  
  if (!is.null(adaptor)) adaptor = validate_function_or_adaptor(obj = adaptor)
  
  
  
  #-------------------------------------
  # First, if both are NULL => no replication
  #-------------------------------------
  
  #  Combine each method: K, K1, K2, etc.
  #  We do a simple for loop accumulation.
  #  For scalar results (like K, K3operator, etc.), we sum them.
  #  For vector results (K1, etc.), we do elementwise sum.
  #  For matrix results (K2, etc.), we do matrix sum.
  #  For constraints, we concatenate them, etc.
  
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
      print(cg$K4operatorAABB(tvec, param, Q1, Q2))
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
  
  res <- createCGF(
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
  
  if (!is.null(adaptor)) res <- adaptCGF(cgf = res, param_adaptor = adaptor)
  if (is.null(block_size) && is.null(iidReps)) return(res)
  if (!is.null(iidReps) && iidReps == 1) return(res)
  
  iidReplicatesCGF(cgf = res, iidReps = iidReps, block_size = block_size, ...)
}
