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

  child_K2_factor <- .K2_factor_method(cgf)

  # ----------------------------------------------------------------
  #   Wrap the five required CGF methods
  # ----------------------------------------------------------------

  K <- function(tvec, param) {
    cgf$K(tvec, param_adaptor(param))
  }
  K1 <- function(tvec, param) {
    cgf$K1(tvec, param_adaptor(param))
  }
  K2 <- function(tvec, param) {
    cgf$K2(tvec, param_adaptor(param))
  }
  K3operator <- function(tvec, param, v1, v2, v3) {
    cgf$K3operator(tvec, param_adaptor(param), v1, v2, v3)
  }
  K4operator <- function(tvec, param, v1, v2, v3, v4) {
    cgf$K4operator(tvec, param_adaptor(param), v1, v2, v3, v4)
  }

  # ----------------------------------------------------------------
  #    Wrap optional public methods
  #    (like K2operator, K2operatorAK2AT, K4operatorAABB, etc.)
  # ----------------------------------------------------------------
  K2operator <- function(tvec, param, x, y) {
    cgf$K2operator(tvec, param_adaptor(param), x, y)
  }


  K2operatorAK2AT <- function(tvec, param, A) {
    cgf$K2operatorAK2AT(tvec, param_adaptor(param), A)
  }

  K2_factor <- NULL
  if (!is.null(child_K2_factor)) {
    K2_factor <- function(tvec, param, A) {
      child_K2_factor(tvec, param_adaptor(param), A)
    }
  }

  K4operatorAABB <- function(tvec, param, Q) {
    cgf$K4operatorAABB(tvec, param_adaptor(param), Q)
  }

  K3K3operatorAABBCC <- function(tvec, param, Q) {
    cgf$K3K3operatorAABBCC(tvec, param_adaptor(param), Q)
  }

  K3K3operatorABCABC <- function(tvec, param, Q) {
    cgf$K3K3operatorABCABC(tvec, param_adaptor(param), Q)
  }

  ineq_constraint <- function(tvec, param) {
    cgf$ineq_constraint(tvec, param_adaptor(param))
  }

  analytic_tvec_hat <- NULL
  if (isTRUE(cgf$has_analytic_tvec_hat)) {
    child_analytic_tvec_hat <- cgf$.private_api$analytic_tvec_hat_func
    analytic_tvec_hat <- function(x, param) {
      child_analytic_tvec_hat(x, param_adaptor(param))
    }
  }

  base_tilting_exponent <- cgf$.private_api$tilting_exponent
  base_neg_ll <- cgf$.private_api$neg_ll
  base_func_T <- cgf$.private_api$func_T
  base_K4operatorAABB_factored     <- cgf$.private_api$K4operatorAABB_factored
  base_K3K3operatorAABBCC_factored <- cgf$.private_api$K3K3operatorAABBCC_factored
  base_K3K3operatorABCABC_factored <- cgf$.private_api$K3K3operatorABCABC_factored

  tilting_exponent <- function(tvec, param) base_tilting_exponent(tvec, param_adaptor(param))
  neg_ll <- function(tvec, param) base_neg_ll(tvec, param_adaptor(param))
  func_T <- function(tvec, param) base_func_T(tvec, param_adaptor(param))
  K4operatorAABB_factored     <- function(tvec, param, A, d) base_K4operatorAABB_factored(tvec, param_adaptor(param), A, d)
  K3K3operatorAABBCC_factored <- function(tvec, param, A, d) base_K3K3operatorAABBCC_factored(tvec, param_adaptor(param), A, d)
  K3K3operatorABCABC_factored <- function(tvec, param, A, d) base_K3K3operatorABCABC_factored(tvec, param_adaptor(param), A, d)
  K4operatorAABB_factored <- .factored_delegate_mark(
    K4operatorAABB_factored,
    .factored_delegate_is_safe(base_K4operatorAABB_factored)
  )
  K3K3operatorAABBCC_factored <- .factored_delegate_mark(
    K3K3operatorAABBCC_factored,
    .factored_delegate_is_safe(base_K3K3operatorAABBCC_factored)
  )
  K3K3operatorABCABC_factored <- .factored_delegate_mark(
    K3K3operatorABCABC_factored,
    .factored_delegate_is_safe(base_K3K3operatorABCABC_factored)
  )







  K2_solve <- function(tvec, param, rhs) {
    cgf$K2_solve(tvec, param_adaptor(param), rhs)
  }

  logdetK2 <- function(tvec, param) {
    cgf$logdetK2(tvec, param_adaptor(param))
  }

  rsim <- NULL
  if (isTRUE(cgf$has_rsim)) {
    # This becomes *rsim_func* in the wrapped CGF.
    # Always return a matrix; flattening is handled by the public CGF$rsim().
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      cgf$rsim(
        n = n,
        vector_length = vector_length,
        parameter_vector = param_adaptor(parameter_vector),
        tvec = tvec,
        flatten = FALSE,
        ...
      )
    }
  }

  op_name <- paste0("A-{", cgf$call_history, "}")

  structured_pair_safe <- .K2_structured_pair_is_safe(cgf)
  K2_solve <- .K2_structured_pair_mark(K2_solve, structured_pair_safe)
  logdetK2 <- .K2_structured_pair_mark(logdetK2, structured_pair_safe)

  # ----------------------------------------------------------------
  # New CGF using createCGF(), passing these wrappers
  # ----------------------------------------------------------------

  cgf_args <- list(
    K = K,
    K1 = K1,
    K2 = K2,
    K3operator = K3operator,
    K4operator = K4operator,
    tilting_exponent = tilting_exponent,
    neg_ll = neg_ll,
    func_T = func_T,
    ineq_constraint = ineq_constraint,
    analytic_tvec_hat = analytic_tvec_hat,
    K2operator = K2operator,
    K2operatorAK2AT = K2operatorAK2AT,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    K2_solve = K2_solve,
    logdetK2 = logdetK2,
    K2_factor = K2_factor,
    rsim = rsim,
    op_name = op_name
  )

  do.call(createCGF, c(cgf_args, list(...)))
}






#' Create an adapted CGF object using an adaptor function
#'
#' @description
#' Constructs a new \code{CGF} object by adapting the parameter vector using a user-supplied
#' \code{adaptor} before invoking the original \code{CGF}'s methods.
#'
#'
#' @param cgf A `CGF` object to be adapted.
#' @param adaptor An \code{adaptor} or a function with signature \code{function(theta) -> adapted_param}.
#' @param ... Additional named arguments passed to \code{\link{createCGF}}, the `CGF` object creation function (rarely needed).
#'
#'
#' @details
#' This function is useful when you have a \code{CGF} that expects its parameter vector
#' to be in a certain format, but your high-level model provides a different parameter
#' structure. By specifying a \code{adaptor}, you can dynamically
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
#' # adaptor that converts 'theta' to c(lambda1, lambda2)
#' mapThetaToDistParams <- function(theta) c(theta, 2*theta)
#'
#' # Now adapt sum_cgf so it only needs a single 'theta':
#' adapted_cgf <- adaptCGF(cgf = sum_cgf, adaptor = mapThetaToDistParams)
#' theta0 <- 5
#' adapted_cgf$K1(0, theta0)
#' # OR sum_cgf$K1(0, c(theta0, 2*theta0))
#' # OR PoissonCGF$K1(0, 3*theta0)
#' }
#'
#' @export
adaptCGF <- function(cgf, adaptor, ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be a CGF object (inherits from 'CGF').")
  param_adaptor <- validate_function_or_adaptor(adaptor)
  .adaptCGF_internal(cgf, param_adaptor, ...)
}
