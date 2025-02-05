# R/D-SumsOfiidCGF.R


.sumOfiidCGF_internal <- function(cgf,
                                  n_fn,
                                  ...) {
  
  
  # K(tvec,param) = n * old_cgf$K(tvec,param)
  Kfun <- function(tvec, param) n_fn(param)*cgf$K(tvec, param)
  K1fun <- function(tvec, param) n_fn(param)*cgf$K1(tvec, param)
  K2fun <- function(tvec, param) {
    ###### issues when "cgf$K2(tvec, param)" is 1 by 1 matrix
    ###### 4*RTMB::advector(matrix(3)) returns a vector and not a matrix
    ###### for now we handle this separately
    n <- n_fn(param)
    k2_val <- cgf$K2(tvec, param)
    if(length(k2_val) == 1) return(n %*% k2_val)
    n*cgf$K2(tvec, param)
  }
  K3opfun <- function(tvec, param, v1, v2, v3) n_fn(param)*cgf$K3operator(tvec, param, v1, v2, v3)
  K4opfun <- function(tvec, param, v1, v2, v3, v4) n_fn(param)*cgf$K4operator(tvec, param, v1, v2, v3, v4)
  
  
  #    tilting_exponent(tvec,param) = n * old_cgf$tilting_exponent(tvec,param)
  #    by definition. 
  #    tilting_exponent = K - sum(t*K1). The factor n naturally emerges from each part.
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  tiltingfun <- function(tvec, param) n_fn(param)*tilting_exponent(tvec, param)
  
  
  
  
  
  K2opfun <- function(tvec, param, x, y) n_fn(param)*cgf$K2operator(tvec, param, x, y)
  K2opAK2ATfun <- function(tvec, param, A) n_fn(param)*cgf$K2operatorAK2AT(tvec,param, A)
  
  
  #   func_T(tvec,param) => base cgf's T / n
  funcT <- cgf$.get_private_method("func_T")
  func_Tfun <- function(tvec, param) funcT(tvec, param)/n_fn(param)
  
  K3K3operatorAABBCC_factored <- cgf$.get_private_method("K3K3operatorAABBCC_factored")
  K4operatorAABB_factored <- cgf$.get_private_method("K4operatorAABB_factored")
  K3K3operatorABCABC_factored <- cgf$.get_private_method("K3K3operatorABCABC_factored")
  
  K3K3AABBCCfactored_fun <- function(tvec, param, A1, d1, A2, d2, A3, d3) {
    n <- n_fn(param)
    K3K3operatorAABBCC_factored(tvec, param, 
                                A1, n*d1,
                                A2, n*d2,
                                A3, n*d3)/n
  }
  K4AABBfactored_fun <- function(tvec, param, A1, d1, A2, d2) {
    n <- n_fn(param)
    K4operatorAABB_factored(tvec, param, 
                            A1, n*d1,
                            A2, n*d2)/n
  }
  K3K3ABCABCfactored_fun <- function(tvec, param, A1, d1, A2, d2, A3, d3) {
    n <- n_fn(param)
    K3K3operatorABCABC_factored(tvec, param, 
                                A1, n*d1,
                                A2, n*d2,
                                A3, n*d3)/n
  }
  
  K3K3AABBCCfun <- function(tvec, param, Q1, Q2, Q3) n_fn(param)^2*cgf$K3K3operatorAABBCC(tvec, param, Q1, Q2, Q3)
  K3K3ABCABCfun <- function(tvec, param, Q1, Q2, Q3) n_fn(param)^2*cgf$K3K3operatorABCABC(tvec, param, Q1, Q2, Q3)
  K4AABBfun <- function(tvec, param, Q1, Q2) n_fn(param)*cgf$K4operatorAABB(tvec, param, Q1, Q2)
  
  
  
  ineqfun <- function(tvec, param) cgf$ineq_constraint(tvec, param)
  
  
  # for the analytic_tvec_hat_func:
  # new_analytic(y, param) = old_analytic( y/n, param )
  if ( cgf$has_analytic_tvec_hat() ) {
    analytic_tvec_hat_func <- function(y, param) {
      cgf$analytic_tvec_hat(y/n_fn(param), param)
    }
  } else {
    # If the base CGF had no valid function, just return NULL
    analytic_tvec_hat_func <- NULL
  }
  
  
  
  
  #     Build the new CGF using createCGF
  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    
    tilting_exponent = tiltingfun,
    # neg_ll           = negllfun, ##### We use the default neg_ll, but we could define a new one???.
    
    K2operator       = K2opfun,
    K2operatorAK2AT  = K2opAK2ATfun,
    K4operatorAABB   = K4AABBfun,
    K3K3operatorAABBCC = K3K3AABBCCfun,
    K3K3operatorABCABC = K3K3ABCABCfun,
    K4operatorAABB_factored     = K4AABBfactored_fun,
    K3K3operatorAABBCC_factored = K3K3AABBCCfactored_fun,
    K3K3operatorABCABC_factored = K3K3ABCABCfactored_fun,
    func_T           = func_Tfun,
    
    ineq_constraint  = ineqfun,
    
    analytic_tvec_hat_func = analytic_tvec_hat_func,
    
    op_name       = c(cgf$call_history, "sumOfiidCGF"),
    ...
  )
}







#' Create a CGF object for the sum of `n` i.i.d. random variables.
#'
#' @description
#' Given a base CGF object `cgf` describing a random variable \eqn{X},
#' this function returns a new CGF object for the sum of `n` i.i.d. copies of \eqn{X},
#' i.e. `Y = X1 + X2 + ... + Xn`.
#' 
#' The argument \code{n} can be:
#' \itemize{
#'   \item A positive numeric scalar (e.g., \code{n=5}).
#'   \item A function or adaptor, \eqn{\theta \mapsto n(\theta)},
#'         returning a positive numeric value for some parameter vector \code{theta}.
#' }
#' 
#' @param cgf A `CGF` object describing the base random variable \eqn{X}.
#' @param n Either a positive scalar or a function(\code{theta}) -> (positive scalar).
#' @param block_size Either \code{NULL} or a positive integer specifying the block size
#'   for i.i.d. replication of the resulting CGF. Defaults to \code{NULL}.
#' @param iidReps Either \code{NULL} or a positive integer specifying how many i.i.d. blocks will be expected by the resulting CGF. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to the `CGF` creation function
#' 
#'
#' @return  A new `CGF` object for \eqn{\;Y = \sum_{i=1}^{n(\theta)} X_i.}
#' @examples
#' \dontrun{
#' # Example: If 'cgf' describes a single Poisson(\lambda),
#' # then sumOfiidCGF(cgf, n=5) describes Poisson(5*lambda).
#'
#' # Evaluate K, K1, etc. on the new object at t=0.1
#' # new_cgf <- sumOfiidCGF(cgf, n=5)
#' # new_cgf$K(0.1, param= c(...)) # => 5 times the old K(0.1)
#' }
#' @export
sumOfiidCGF <- function(cgf, n, 
                        block_size = NULL,
                        iidReps    = NULL,
                        ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be an object of class 'CGF'.")
  
  if (is.numeric(n)) {
    if (length(n) != 1 || n <= 0) {
      stop("'n' must be a positive numeric scalar.")
    }
    n_fn <- function(theta) n+0*theta[1] # a trivial function
  } else {
    # Otherwise interpret n as an adaptor/function
    n_fn <- validate_function_or_adaptor(n)
  }
  
  base_cgf <- .sumOfiidCGF_internal(cgf, n_fn, ...)
  
  
  if (is.null(block_size) && is.null(iidReps)) return(base_cgf)
  if (!is.null(iidReps) && iidReps == 1) return(base_cgf)
  
  iidReplicatesCGF(cgf = base_cgf, iidReps = iidReps, block_size = block_size)
}
