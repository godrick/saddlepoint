# --------------------------------------------------------------------
# File: R/D-AdaptCGF.R
#
# This file provides a factory function `adaptCGF()` that takes:
#   1) An existing CGF object (`cgf`)
#   2) A parameter adaptor function (`param_adaptor`)
# and returns a *new* CGF object whose methods adapt parameters before
# delegating to `cgf`.
#
# Usage:
#   adapted <- adaptCGF(original_cgf, function(param) { param * 2 })
#   # Then calling adapted$K(t, param) uses param*2 internally.
#
# Implementation Notes:
#  - We call `createCGF()` to build a fresh CGF. 
#  - For each method, we define a new function that calls `cgf$Method(tvec, param_adaptor(param), ...)`.
#  - Because we do not modify the original object (nor inherit from CGF),
#    we avoid "locked binding" issues entirely.
# --------------------------------------------------------------------




.adaptCGF_internal <- function(cgf, param_adaptor, ...){
  
  # ----------------------------------------------------------------
  #   Wrap the five required CGF methods
  # ----------------------------------------------------------------
  
  wrapped_K <- function(tvec, param) {
    cgf$K(tvec, param_adaptor(param))
  }
  wrapped_K1 <- function(tvec, param) {
    cgf$K1(tvec, param_adaptor(param))
  }
  wrapped_K2 <- function(tvec, param) {
    cgf$K2(tvec, param_adaptor(param))
  }
  wrapped_K3operator <- function(tvec, param, v1, v2, v3) {
    cgf$K3operator(tvec, param_adaptor(param), v1, v2, v3)
  }
  wrapped_K4operator <- function(tvec, param, v1, v2, v3, v4) {
    cgf$K4operator(tvec, param_adaptor(param), v1, v2, v3, v4)
  }
  
  # ----------------------------------------------------------------
  #    Wrap optional public methods
  #    (like K2operator, K2operatorAK2AT, K4operatorAABB, etc.)
  # ----------------------------------------------------------------
  wrapped_K2operator <- function(tvec, param, x, y) {
    cgf$K2operator(tvec, param_adaptor(param), x, y)
  }
  
  
  wrapped_K2operatorAK2AT <- function(tvec, param, A) {
    cgf$K2operatorAK2AT(tvec, param_adaptor(param), A)
  }
  
  wrapped_K4operatorAABB <- function(tvec, param, Q1, Q2) {
    cgf$K4operatorAABB(tvec, param_adaptor(param), Q1, Q2)
  }
  
  wrapped_K3K3operatorAABBCC <- function(tvec, param, Q1, Q2, Q3) {
    cgf$K3K3operatorAABBCC(tvec, param_adaptor(param), Q1, Q2, Q3)
  }
  
  wrapped_K3K3operatorABCABC <- function(tvec, param, Q1, Q2, Q3) {
    cgf$K3K3operatorABCABC(tvec, param_adaptor(param), Q1, Q2, Q3)
  }
  
  wrapped_ineq_constraint <- function(tvec, param) {
    cgf$ineq_constraint(tvec, param_adaptor(param))
  }
  
  if(cgf$has_analytic_tvec_hat()) {
    wrapped_analytic_tvec_hat_func <- function(x, param) {
      cgf$analytic_tvec_hat(x, param_adaptor(param))
    }
  } else {
    wrapped_analytic_tvec_hat_func <- NULL
  }
  
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  neg_ll <- cgf$.get_private_method("neg_ll")
  func_T <- cgf$.get_private_method("func_T")
  K4operatorAABB_factored     <- cgf$.get_private_method("K4operatorAABB_factored")
  K3K3operatorAABBCC_factored <- cgf$.get_private_method("K3K3operatorAABBCC_factored")
  K3K3operatorABCABC_factored <- cgf$.get_private_method("K3K3operatorABCABC_factored")
  
  wrapped_tilting_exponent <- function(tvec, param) tilting_exponent(tvec, param_adaptor(param))
  wrapped_neg_ll <- function(tvec, param) neg_ll(tvec, param_adaptor(param))
  wrapped_func_T <- function(tvec, param) func_T(tvec, param_adaptor(param))
  wrapped_K4operatorAABB_factored     <- function(tvec, param, A1, d1, A2, d2) K4operatorAABB_factored(tvec, param_adaptor(param), A1, d1, A2, d2)
  wrapped_K3K3operatorAABBCC_factored <- function(tvec, param, A1, d1, A2, d2, A3, d3) K3K3operatorAABBCC_factored(tvec, param_adaptor(param), A1, d1, A2, d2, A3, d3)
  wrapped_K3K3operatorABCABC_factored <- function(tvec, param, A1, d1, A2, d2, A3, d3) K3K3operatorABCABC_factored(tvec, param_adaptor(param), A1, d1, A2, d2, A3, d3)
  
  op_name = paste0("A-{", cgf$call_history, "}")
  
  # ----------------------------------------------------------------
  # New CGF using createCGF(), passing these wrappers
  # ----------------------------------------------------------------
  
  createCGF(
    K = wrapped_K,
    K1 = wrapped_K1,
    K2 = wrapped_K2,
    K3operator = wrapped_K3operator,
    K4operator = wrapped_K4operator,
    ineq_constraint = wrapped_ineq_constraint,
    analytic_tvec_hat_func = wrapped_analytic_tvec_hat_func,
    op_name = op_name,   
    
    tilting_exponent = wrapped_tilting_exponent,
    neg_ll = wrapped_neg_ll,
    func_T = wrapped_func_T,
    K4operatorAABB = wrapped_K4operatorAABB,
    K3K3operatorAABBCC = wrapped_K3K3operatorAABBCC,
    K3K3operatorABCABC = wrapped_K3K3operatorABCABC,
    K4operatorAABB_factored = wrapped_K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = wrapped_K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = wrapped_K3K3operatorABCABC_factored,
    K2operator = wrapped_K2operator,
    K2operatorAK2AT = wrapped_K2operatorAK2AT,
    
    ...
  )
}






#' Create an adapted CGF object using an adaptor function
#'
#' @description
#' Constructs a new \code{CGF} object by adapting the parameter vector using a user-supplied
#' \code{param_adaptor} before invoking the original \code{CGF}'s methods.
#'
#'
#' @param cgf A `CGF` object to be adapted.
#' @param param_adaptor An \code{adaptor} or a function with signature \code{function(theta) -> adapted_param}.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, the `CGF` object creation function (rarely needed).
#' 
#' 
#' @details
#' This function is useful when you have a \code{CGF} that expects its parameter vector
#' to be in a certain format, but your high-level model provides a different parameter
#' structure. By specifying a \code{param_adaptor}, you can dynamically
#' translate your model parameters to those that \code{cgf} requires.
#'
#' @return A new `CGF` object.
#'
#' @examples
#' \dontrun{
#' ## Example: Suppose you have a sum of two Poisson r.v.s, each with
#' ##   a different lambda, but your model uses a single parameter 'theta'
#' ##   from which the two lambdas are derived (lambda1 = theta, lambda2 = 2*theta).
#'
#' ## Scenario:
#' ## Y1 ~ Poisson(theta)
#' ## Y2 ~ Poisson(2*theta)
#' ## Y = Y1 + Y2 ~ Poisson(3*theta)
#' ## Model Parameters: theta
#' ## CGF for Y expects: distribution_params = c(lambda1, lambda2) = c(theta, 2*theta)
#' 
#'
#' # First, the individual Poisson CGFs with separate lambda parameters
#' # These will expect a vector of length 2 for the two lambdas
#' K_Y1 <- PoissonModelCGF(lambda = adaptor(indices = 1))   # For Y1: lambda1 => the first parameter
#' K_Y2 <- PoissonModelCGF(lambda = adaptor(indices = 2))   # For Y2: lambda2 => the second parameter
#' 
#' # sum_cgf is a CGF that expects 2 parameters for the two Poisson variables
#' sum_cgf <- sumOfIndependentCGF(cgf_list = list(K_Y1, K_Y2))
#' 
#' # param_adaptor that converts 'theta' to c(lambda1, lambda2)
#' mapThetaToDistParams <- function(theta) c(theta, 2*theta) 
#' 
#' # Now adapt sum_cgf so it only needs a single 'theta':
#' adapted_cgf <- adaptCGF(cgf = sum_cgf, param_adaptor = mapThetaToDistParams)
#' theta0 <- 5
#' adapted_cgf$K1(0, theta0)  
#' # OR sum_cgf$K1(0, c(theta0, 2*theta0))
#' # OR PoissonCGF$K1(0, 3*theta0)
#' }
#'
#' @export
adaptCGF <- function(cgf, param_adaptor, ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object (inherits from 'CGF').")
  param_adaptor <- validate_function_or_adaptor(param_adaptor)
  .adaptCGF_internal(cgf, param_adaptor, ...)
}
