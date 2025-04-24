# R/NormalCGF.R
# Objects: NormalCGF, MultivariateNormalModelCGF






#' Univariate Normal CGF object
#'
#' A ready-to-use CGF object for the univariate Normal distribution with mean
#' \eqn{\mu} and standard deviation \eqn{\sigma}. By default, \code{NormalCGF} evaluation for i.i.d. replicates.
#' If \code{tvec} has length \eqn{n}, we interpret that as \eqn{n} i.i.d. replicates. 
#'
#' @details
#' **Parameter Vector**: The \code{parameter_vector} used when calling methods such as `K(tvec, parameter_vector)`
#' must be a numeric vector of the form \eqn{(\mu, \sigma)} in that order.
#' 
#'
#' @format
#' An object of class \code{CGF} (R6), with the usual methods \code{K, K1, K2, etc.}
#'
#' @examples
#' NormalCGF$K1(0, c(5, 0.5)) # (expected value for N(5, 0.5) is 5)
#'
#' @export
NormalCGF <- createCGF_fromVectorisedFunctions(
  
  # K(tvec) = sum_k [ mu t_k + 0.5 sigma^2 t_k^2 ]
  K_vectorized_func = function(tvec, param) {
    mu    <- param[1]
    sigma2 <- param[2]^2
    mu*tvec + 0.5*sigma2*tvec^2
  },
  # K1(tvec) => for each component i: mu + sigma^2 t_i
  K1_vectorized_func = function(tvec, param) {
    mu    <- param[1]
    sigma2 <- param[2]^2
    mu + sigma2 * tvec
  },
  # K2(tvec) => a diagonal matrix with diagonal = sigma^2
  # For a length-n tvec, we expect an n-by-n matrix = sigma^2*I (in the parent class).
  K2_vectorized_func = function(tvec, param) {
    sigma2 <- param[2]^2
    rep(sigma2, length(tvec))
  },
  # K3 is identically 0, but we must define a "vectorized" version that 
  # takes (tvec, param) and acts on w1, w2, w3 => always 0
  K3_vectorized_func = function(tvec, param) {
    # 0 for any input => sum(0)
    rep(0, length(tvec))  
  },
  # K4 is identically 0 as well => 0
  K4_vectorized_func = function(tvec, param) {
    rep(0, length(tvec))
  },
  analytic_tvec_hat_func = function(y, param) { (y - param[1])/(param[2]^2)  },
  op_name = "NormalCGF"
)










.MultivariateNormalModelCGF_internal <- function(iidReps, ...) {
    Kfun = function(tvec, param) {
      d <- (-1 + sqrt(1 + 4*length(param))) / 2
      if (length(tvec) != d*iidReps) stop(sprintf("tvec has length %d, but expecting d*iidReps = %d*%d", length(tvec), d, iidReps))
      muVec     <- param[1:d]
      Sigma_mat <- matrix(param[(d+1):length(param)], nrow=d, ncol=d, byrow=FALSE)  ##### Matrix(advector(...),...) yields an error
      
      tmat <- matrix(tvec, nrow = d)
      block_vals <- as.vector(crossprod(muVec, tmat)) + 0.5 * colSums(tmat * (Sigma_mat %*% tmat))
      sum(block_vals)
    }
    
    K1fun <- function(tvec, param) {
      d <- (-1 + sqrt(1 + 4*length(param))) / 2
      if (length(tvec) != d*iidReps) stop(sprintf("tvec has length %d, but expecting d*iidReps = %d*%d", length(tvec), d, iidReps))
      muVec     <- param[1:d]
      Sigma_mat <- matrix(param[(d+1):length(param)], nrow=d, ncol=d, byrow=FALSE)
      
      tmat <- matrix(tvec, nrow = d)
      as.vector(muVec + (Sigma_mat %*% tmat))
    }
    
    K2fun <- function(tvec, param) {
      d <- (-1 + sqrt(1 + 4*length(param))) / 2
      if (length(tvec) != d*iidReps) stop(sprintf("tvec has length %d, but expecting d*iidReps = %d*%d", length(tvec), d, iidReps))
      Sigma_mat <- matrix(param[(d+1):length(param)], nrow=d, ncol=d, byrow=FALSE)
      big_mat   <- Matrix::Matrix(0, nrow = length(tvec), ncol = length(tvec)) * param[1] #### decide on sparse?
      for (i in seq_len(iidReps)) {
        idx <- seq.int((i - 1)*d + 1, i*d)
        big_mat[idx, idx] <- Sigma_mat
      }
      big_mat
    }
    
    K3opfun <- function(tvec, param, v1, v2, v3) 0
    K4opfun <- function(tvec, param, v1, v2, v3, v4) 0
    
    
    K4AABBfun <- function(tvec, param, Q1, Q2) 0
    K3K3operatorAABBCCfun <- function(tvec, param, Q1, Q2, Q3) 0
    K3K3operatorABCABCfun <- function(tvec, param, Q1, Q2, Q3) 0
    
    
    K3K3operatorABCABC_factoredfun <- function(tvec, param, A1, d1, A2, d2, A3, d3) 0
    K3K3operatorAABBCC_factoredfun <- function(tvec, param, A1, d1, A2, d2, A3, d3) 0
    K4operatorAABB_factoredfun     <- function(tvec, param, A1, d1, A2, d2) 0
    
    func_Tfun <- function(tvec, param) 0
  
    saddlepoint_t_MVN <- function(y, param) {
      d <- (-1 + sqrt(1 + 4*length(param))) / 2
      if (length(y) != d*iidReps) stop(sprintf("y has length %d, but expecting d*iidReps = %d*%d", length(y), d, iidReps))
      muVec     <- param[1:d]
      Sigma_mat <- matrix(param[(d+1):length(param)], nrow=d, ncol=d, byrow=FALSE)
      y_mat <- matrix(y, nrow = d)
      res_mat <- solve(Sigma_mat, y_mat - muVec)
      as.vector(res_mat)
    }
    
    createCGF(
      K = Kfun, 
      K1 = K1fun, 
      K2 = K2fun, 
      K3operator = K3opfun, 
      K4operator = K4opfun, 
      analytic_tvec_hat_func = saddlepoint_t_MVN,
      op_name = "MultivariateNormalModelCGF",
      
      func_T = func_Tfun,
      K4operatorAABB = K4AABBfun,
      K3K3operatorAABBCC = K3K3operatorAABBCCfun,
      K3K3operatorABCABC = K3K3operatorABCABCfun,
      K4operatorAABB_factored = K4operatorAABB_factoredfun,
      K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factoredfun,
      K3K3operatorABCABC_factored = K3K3operatorABCABC_factoredfun,
      ...
    )
    
}





#' Create a Multivariate Normal CGF Object 
#'
#' @description
#' Creates a CGF for an arbitrary \code{d}-dimensional normal distribution.
#' The user supplies two functions:
#' \itemize{
#'   \item \code{mu(theta)} returning a length-\eqn{d} numeric vector.
#'   \item \code{sigma(theta)} returning a \eqn{d \times d} matrix.
#' }
#' 
#' @details
#' **I.I.D. Replicates**:
#' By setting \code{iidReps} to a positive integer \eqn{m}, you declare that the input vector \eqn{tvec} will be split
#' into \eqn{m} blocks of equal size, each corresponding to one i.i.d. multivariate normal sample. 
#' If \code{iidReps} is \code{"any"}, no length restriction is enforced on \eqn{tvec}, allowing flexible usage.
#'
#' @param mu A function(\code{theta}) -> numeric vector (the mean).
#' @param sigma A function(\code{theta}) -> matrix (the covariance).
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional arguments passed to \code{\link{createCGF}}
#'
#'
#'
#' @return A `CGF` object.
#'
#' @export
MultivariateNormalModelCGF <- function(mu, sigma, iidReps = "any", ...) {
  
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  
  mu_fn <- validate_function_or_adaptor(mu)
  sigma_fn <- validate_function_or_adaptor(sigma)
  
  # param_adaptor: calls mu(theta) and sigma(theta), flattens
  param_adaptor_ <- function(theta) {
    muVal    <- mu_fn(theta)      # must be length d
    sigmaVal <- sigma_fn(theta)   # must be d x d
    if (nrow(sigmaVal) != length(muVal) || ncol(sigmaVal) != length(muVal)) stop("sigma(theta) must be a square matrix")
    c( muVal, as.vector(sigmaVal) )  # sigmaVal is column-major => must unwrap consistently (byrow = FALSE)
  }
  
  base_cgf <- .MultivariateNormalModelCGF_internal(iidReps, ...)
  adaptCGF(cgf = base_cgf, adaptor = param_adaptor_)
}






