# R/NormalCGF.R
# Objects: NormalCGF, UnivariateNormalModelCGF, MultivariateNormalModelCGF

#' Univariate Normal CGF object
#'
#' A ready-to-use CGF object for the univariate Normal distribution with mean
#' \eqn{\mu} and standard deviation \eqn{\sigma}. 
#'
#' @details
#' **Parameter Vector**: The \code{parameter_vector} should be 
#' \eqn{(\mu, \sigma)} in that order. 
#' 
#' If \code{tvec} has length \eqn{n}, we interpret that as 
#' \eqn{n} i.i.d. replicates. 
#' 
#'
#' @format
#' An object of class \code{CGF} (R6), with the usual methods \code{K, K1, K2, etc.}
#'
#' @examples
#' # Example: Evaluate K at tvec=c(0.1, -0.2) for mu=2, var=4
#' # param = c(2, 4)
#' #### TODO: Add example
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
  param_adaptor = function(x) x[1:2],
  
  op_name = "NormalCGF"
)




#' Create a Parametric Gaussian CGF Object
#'
#' @description
#' Constructs a CGF object for a univariate normal distribution where mean \eqn{\mu(\theta)} 
#' and variance \eqn{\sigma(\theta)} are derived from user-supplied logic.
#' The user supplies a function that maps their full parameter vector 
#' \code{theta} to \eqn{(\mu, \sigma)}.
#'
#' @param mean A function that accepts a single parameter vector \code{theta} and returns \eqn{\mu}. 
#' @param sd A function that accepts a single parameter vector \code{theta} and returns \eqn{\sigma}. 
#' @param iidReps A positive integer. If \code{iidReps>1}, 
#'   we expect \code{tvec} to have that length.
#' @param ... Additional arguments passed to base CGF creation function.
#'
#'
#' @return A CGF object.
#' @export
UnivariateNormalModelCGF <- function(mean, sd, iidReps = 1, ...) {
  if (!is.function(mean) || !is.function(sd)) stop("`mean` and `sd` must be functions (theta) -> numeric")
  if (!is.numeric(iidReps) || length(iidReps)!=1 || iidReps<1) stop("`iidReps` must be a positive integer.")
  if (length(formals(mean))!=1 || length(formals(sd))!=1) stop("`mean` and `sd` must be functions of a single parameter vector.")
  
  
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) {
      stop(sprintf("`tvec` must have length %d (got %d).", iidReps, length(tvec)))
    }
  }
  
  
  createCGF_fromVectorisedFunctions(
    # K(tvec)
    K_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      param[1]*tvec + 0.5*(param[2]^2)*tvec^2
    },
    # K1(tvec)
    K1_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      mu    <- param[1]
      sigma2 <- param[2]^2
      mu + sigma2 * tvec
    },
    # K2(tvec)
    K2_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      sigma2 <- param[2]^2
      rep(sigma2, length(tvec))
    },
    # K3(tvec)=0
    K3_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      # Return 0 => either a single 0 or a length=1 numeric, but for the operator,
      # we might expect a scalar. We'll just do 0:
      0
    },
    # K4(tvec)=0
    K4_vectorized_func = function(tvec, param) {
      check_tvec(tvec)
      0
    },
    # analytic_tvec_hat: Solve mu + sigma^2 t = y => t= (y - mu)/ sigma^2
    analytic_tvec_hat_func = function(y, param) {
      mu    <- param[1]
      sigma2 <- param[2]^2
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      (y - mu)/sigma2
    },
    param_adaptor = function(theta) c(mean(theta), sd(theta)),
    op_name = "UnivariateNormalModelCGF",
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
#' @param mu A function(\code{theta}) -> numeric vector (the mean).
#' @param sigma A function(\code{theta}) -> matrix (the covariance).
#' @param iidReps An integer. If >1, we replicate blocks for i.i.d. usage.
#' @param ... Additional arguments passed to \code{\link{createCGF}} or 
#'   to \code{\link{iidReplicatesCGF}} if \code{iidReps>1}.
#'
#' @details
#' **Parameter Flattening**: 
#' If \code{mu(theta)} is length \eqn{d}, and \code{Sigma(theta)} is \eqn{d x d},
#' then we store them as \code{c(mu, as.vector(Sigma))} in row-major or column-major form 
#' (we just must be consistent in unwrapping).
#'
#'
#' @return A CGF object.
#'
#' @export
MultivariateNormalModelCGF <- function(mu, sigma, iidReps=1, ...) {
  
  if (!is.function(mu) || !is.function(sigma)) stop("`mu` and `sigma` must be functions.")
  if (length(formals(mu))!=1 || length(formals(sigma))!=1) stop("`mu` and `sigma` must be functions of a single parameter vector.")
  if (!is.numeric(iidReps) || length(iidReps)!=1 || iidReps<1) stop("`iidReps` must be a positive integer.")
  
  
  
  # param_adaptor: calls mu(theta) and sigma(theta), flattens
  param_adaptor_ <- function(theta) {
    muVal    <- mu(theta)      # must be length d
    sigmaVal <- sigma(theta)   # must be d x d
    # We won't know d until we see these results or see tvec, so we won't check here
    # but let's flatten them anyway
    c( muVal, as.vector(sigmaVal) )  # assume sigmaVal is row-major => must unwrap consistently
  }
  
  
  Kfun <- function(tvec, param) {
    # unflatten
    up <- unwrapParam(param, length(tvec))
    muVal    <- up$mu
    SigmaVal <- up$sigma
    d <- up$d
    if (length(tvec) != d*iidReps) {
      stop(sprintf("tvec has length %d, but expecting d*iidReps = %d*%d", length(tvec), d, iidReps))
    }
    # For a single block perspective => if iidReps=1 => length(tvec)=d => standard formula
    # But we are a single-block definition => so we expect tvec length= d if iidReps=1,
    # or we let iidReplicatesCGF handle chunking if >1. So we'll interpret this as one block:
    # K(t)= t^T mu + 0.5 t^T sigma t
    # but if length(tvec)= d*iidReps, we still treat it as one big block => meaning a dimension of d*iidReps??? 
    # Actually, let's keep it simpler: if we are the single-block cgf, we want length(tvec)= d, ignoring iidReps. 
    # All checks for iidReps will be done in iidReplicatesCGF.
    as.vector(t(tvec) %*% muVal + 0.5*( t(tvec) %*% SigmaVal %*% tvec ))
  }
  
  K1fun <- function(tvec, param) {
    up <- unwrapParam(param, length(tvec))
    muVal    <- up$mu
    SigmaVal <- up$sigma
    d <- up$d
    if (length(tvec) != d) stop("length of tvec is inconsistent with the parameter dimension")
    as.vector(muVal + (SigmaVal %*% tvec))
  }
  
  K2fun <- function(tvec, param) {
    up <- unwrapParam(param, length(tvec))
    d <- up$d
    if (length(tvec) != d) stop("length of tvec is inconsistent with the parameter dimension")
    up$sigma
  }
  
  K3opfun <- function(tvec, param, v1, v2, v3) 0
  K4opfun <- function(tvec, param, v1, v2, v3, v4) 0
  
  
  K2opfun <- function(tvec, param, x, y) {
    up <- unwrapParam(param, length(tvec))
    as.vector( t(x) %*% (up$sigma %*% y) )
  }
   
  K2opAK2ATfun <- function(tvec, param, A) {
    up <- unwrapParam(param, length(tvec))
    A %*% up$sigma %*% t(A)
  }

  K4AABBfun <- function(tvec, param, Q1, Q2) 0
  K3K3operatorAABBCCfun <- function(tvec, param, Q1, Q2, Q3) 0
  K3K3operatorABCABCfun <- function(tvec, param, Q1, Q2, Q3) 0
  
  
  K3K3operatorABCABC_factoredfun <- function(tvec, param, A1, d1, A2, d2, A3, d3) 0
  K3K3operatorAABBCC_factoredfun <- function(tvec, param, A1, d1, A2, d2, A3, d3) 0
  K4operatorAABB_factoredfun     <- function(tvec, param, A1, d1, A2, d2) 0
  
  func_Tfun <- function(tvec, param) 0
  
  saddlepoint_t_MVN <- function(y, param) {
    up <- unwrapParam(param, length(y))
    solve(up$sigma, y - up$mu)
  }
  
  
  singleBlock_cgf <- createCGF(
                  K = Kfun, 
                  K1 = K1fun, 
                  K2 = K2fun, 
                  K3operator = K3opfun, 
                  K4operator = K4opfun, 
                  param_adaptor = param_adaptor_,
                  analytic_tvec_hat_func = saddlepoint_t_MVN,
                  op_name = "MultivariateNormalModelCGF",
              
                  func_T = func_Tfun,
                  K4operatorAABB = K4AABBfun,
                  K3K3operatorAABBCC = K3K3operatorAABBCCfun,
                  K3K3operatorABCABC = K3K3operatorABCABCfun,
                  K4operatorAABB_factored = K4operatorAABB_factoredfun,
                  K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factoredfun,
                  K3K3operatorABCABC_factored = K3K3operatorABCABC_factoredfun,
                  K2operator = K2opfun,
                  K2operatorAK2AT = K2opAK2ATfun,
                  ...
                  )
  
  if (iidReps==1) return(singleBlock_cgf)

  iidReplicatesCGF(cgf = singleBlock_cgf, iidReps=iidReps, ...)
}



# unflatten => returns list(mu=..., sigma=...)
# We'll also deduce 'd' from param length or from tvec.
#' @noRd
unwrapParam <- function(param, tvecLength) {
  # If single block => tvecLength = d. If multiple => tvecLength= d*iidReps => d= tvecLength / iidReps
  # So let's define:
  d <- tvecLength 
  # param must have length (d + d^2)
  expectedParamLen <- d + d^2
  if (length(param) != expectedParamLen) {
    stop(sprintf("param has length %d; expected %d = d + d*d with d=%.0f based on tvec length %d / 1 %d",
                 length(param), expectedParamLen, d, tvecLength, iidReps))
  }
  muVec  <- param[1:d]
  sigmaF <- param[(d+1) : (d + d*d)]
  # We'll interpret sigmaVal as a matrix(d,d) in row-major form
  sigmaMat <- Matrix(sigmaF, nrow=d, ncol=d, byrow=TRUE)
  
  list(d=d, mu=muVec, sigma=sigmaMat)
}


