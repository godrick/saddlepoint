# R/PoissonCGF.R
# Objects: PoissonCGF, PoissonModelCGF, PoissonNonIdenticalModelCGF





#' Poisson CGF Object
#'
#' A ready-to-use CGF object for the Poisson distribution. 
#' 
#' #' @details
#' **Parameter Vector**: The `parameter_vector` for \code{PoissonCGF} is 
#' interpreted as the rate \eqn{\lambda}. By default, this object is vectorized 
#' for i.i.d. replicates of Poisson variables.
#' 
#' @format An object of class \code{CGF} (an R6 class), with standard methods 
#' \code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.
#'
#' @examples
#' # Evaluate K at t = 0.1 for lambda = 2
#' # PoissonCGF$K(0.1, 2)
#' 
#' @export
PoissonCGF <- createCGF_fromVectorisedFunctions(
  K_vectorized_func = function(tvec, lambda) { lambda * (exp(tvec) - 1) },
  K1_vectorized_func = function(tvec, lambda) {  lambda * exp(tvec) },
  K2_vectorized_func = function(tvec, lambda) { lambda * exp(tvec) },
  K3_vectorized_func = function(tvec, lambda) { lambda * exp(tvec) },
  K4_vectorized_func = function(tvec, lambda) { lambda * exp(tvec) },
  param_adaptor = function(x) x[1]
  ,analytic_tvec_hat_func = function(y, lambda)  log(y / lambda) # solution of the saddlepoint equation
  , op_name = "PoissonCGF"
)


       















#' Create a Parametric Poisson CGF Object
#'
#' @description
#' Creates a CGF object for the Poisson distribution with a rate parameter (\eqn{\lambda}) defined by user-provided parameters.
#' Users supply an adaptor function that maps their parameter vector to \eqn{\lambda}, allowing for flexible parametric modeling.
#'
#' @param lambda A function that accepts a single parameter vector and returns the Poisson rate parameter (\eqn{\lambda}). This function defines how model parameters influence \eqn{\lambda}. 
#' @param iidReps Integer. Specifies the number of i.i.d. replicates expected by the resulting CGF object. Must be a positive integer. Defaults to \code{1}.
#' @param ... Additional arguments passed to the base CGF creation function, such as optional method overrides.
#'
#' @return A CGF object.
#' @export
PoissonModelCGF <- function(lambda, iidReps = 1, ...) {
  # if (!is.function(lambda)) stop("`lambda` must be a valid function that returns the Poisson rate parameter (lambda).")
  #     # lambda should accept only one argument
  # if (length(formals(lambda)) != 1) stop("`lambda` must be a function that accepts exactly one argument.")
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1 || iidReps != as.integer(iidReps)) stop("`iidReps` must be a positive integer.")
  lambda_adaptor <- validate_function_or_adaptor(obj = lambda)
  
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) stop(sprintf("`tvec` must have length %d (got %d).", iidReps, length(tvec)))
  }
  
  # Create the base Poisson CGF
  base_cgf <- createCGF_fromVectorisedFunctions(
    K_vectorized_func  = function(tvec, lam) { 
      check_tvec(tvec)
      lam[1] * (exp(tvec) - 1) 
    },
    K1_vectorized_func = function(tvec, lam) { 
      check_tvec(tvec)
      lam[1] * exp(tvec) 
    },
    K2_vectorized_func = function(tvec, lam) { 
      check_tvec(tvec)
      lam[1] * exp(tvec) 
    },
    K3_vectorized_func = function(tvec, lam) { 
      check_tvec(tvec)
      lam[1] * exp(tvec) 
    },
    K4_vectorized_func = function(tvec, lam) { 
      check_tvec(tvec)
      lam[1] * exp(tvec) 
    },
    param_adaptor          = lambda_adaptor,
    analytic_tvec_hat_func = function(y, lam)  {
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      log(y / lam)
    }, 
    op_name = "PoissonModelCGF",
    ...  # Pass any additional optional arguments
  )
  base_cgf
}











#' A CGF Object for a vector of Non-Identical Poisson variables
#'
#' @description
#' Constructs a CGF object for a \eqn{d}-dimensional vector of independent, non-identically distributed Poisson random variables.
#' Each coordinate \eqn{i} of the vector is associated with its own rate parameter \eqn{\lambda_i}, 
#' determined by the user-provided function \code{lambda_vec}. This function maps a user-specified parameter 
#' vector (e.g., \eqn{\theta}) to distinct rates \eqn{(\lambda_1, \ldots, \lambda_d)}, one for each coordinate.
#'
#' While the primary focus is the vector of independent and non-identical Poisson variables, 
#' the design also supports IID replication of these vectors if needed.
#' 
#'
#' @param lambda_vec A function that accepts a single parameter vector and returns a vector of Poisson rate parameters (\eqn{\lambda_1, \lambda_2, ..., \lambda_d}). This vector defines the rate for each Poisson variable in the vector.
#' @param iidReps Integer. Specifies the number of IID replicate vectors to create. Must be a positive integer.
#' @param ... Additional arguments passed to the base CGF creation function, such as optional method overrides.
#'
#' @return A CGF object.
#' @export
PoissonNonIdenticalModelCGF <- function(lambda_vec, iidReps = 1, ...) {
  
  # if (!is.function(lambda)) stop("`lambda` must be a valid function that returns a vector of Poisson rate parameters.")
  # if (length(formals(lambda)) != 1) stop("`lambda` must be a function that accepts exactly one argument (the parameter vector).")
  lambda_adaptor <- validate_function_or_adaptor(obj = lambda_vec)
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1 || iidReps != as.integer(iidReps)) {
    stop("`iidReps` must be a positive integer.")
  }
  

  base_cgf <- createCGF_fromVectorisedFunctions(
    K_vectorized_func  = function(tvec, params) { 
      ## if (iidReps <= 0 || length(tvec) %% iidReps != 0) stop(sprintf("`tvec` length must be a multiple of `iidReps` (%d) to determine `d`.", iidReps))
      
      # get the number of Poisson variables (d)
      d <- length(params)
      expected_length <- d * iidReps
      # validate tvec length
      if (length(tvec) != expected_length)  stop(sprintf("`tvec` has length %d; expected %d based on parameter length.", length(tvec), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      
      lambdas_replicated * (exp(tvec) - 1) 
    },
    K1_vectorized_func = function(tvec, params) { 
      d <- length(params)
      expected_length <- d * iidReps
      if (length(tvec) != expected_length)  stop(sprintf("`tvec` has length %d; expected %d based on parameter length.", length(tvec), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      lambdas_replicated * exp(tvec) 
    },
    K2_vectorized_func = function(tvec, params) { 
      d <- length(params)
      expected_length <- d * iidReps
      if (length(tvec) != expected_length)  stop(sprintf("`tvec` has length %d; expected %d based on parameter length.", length(tvec), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      lambdas_replicated * exp(tvec)
    },
    K3_vectorized_func = function(tvec, params) { 
      d <- length(params)
      expected_length <- d * iidReps
      if (length(tvec) != expected_length)  stop(sprintf("`tvec` has length %d; expected %d based on parameter length.", length(tvec), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      lambdas_replicated * exp(tvec)
    },
    K4_vectorized_func = function(tvec, params) { 
      d <- length(params)
      expected_length <- d * iidReps
      if (length(tvec) != expected_length)  stop(sprintf("`tvec` has length %d; expected %d based on parameter length.", length(tvec), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      lambdas_replicated * exp(tvec)
    },
    param_adaptor      = lambda_adaptor,
    analytic_tvec_hat_func = function(y, params) {
      d <- length(params)
      expected_length <- d * iidReps
      if (length(y) != expected_length)  stop(sprintf("`y` has length %d; expected %d based on parameter length.", length(y), expected_length))
      lambdas_replicated <- rep(params, times = iidReps)
      log(y / lambdas_replicated)
    },
    op_name = "PoissonNonIdenticalModelCGF",
    ...  # Pass any additional optional arguments
  )
  
  base_cgf
}
