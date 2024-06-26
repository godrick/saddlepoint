#' @title Find maximum likelihood estimates using the saddlepoint likelihood.
#'
#' @description
#' This function uses the nloptr package for optimization to find the maximum likelihood estimates (MLEs). It takes the observed data, a CGF object that represents the CGF of the observed data, and various other parameters related to the optimization procedure.
#'
#' @importFrom nloptr nloptr
#'
#' @param observed.data a numeric vector - See TO DO list
#  TO DO: for a (block.size by m) dataset, nrow(observed.data) represents the
#        dimension of the observed vector and the columns correspond to i.i.d replicates of the vector.
#' @param model.cgf a CGF object corresponding to the observed data.
#' @param starting.theta a numeric vector of starting values for model parameters.
#' @param lb.theta a numeric vector, vector of lower bounds for model parameters, defaults to -Inf for all parameters.
#' @param ub.theta a numeric vector, vector of upper bounds for model parameters, defaults to Inf for all parameters.
#' @param starting.tvec a numeric vector, vector of starting values for saddlepoint parameters, defaults to zeros.
#' @param lb.tvec a numeric vector, vector of lower bounds for saddlepoint parameters, defaults to -Inf for all parameters.
#' @param ub.tvec a numeric vector, vector of upper bounds for saddlepoint parameters, defaults to Inf for all parameters.
#' @param std.error a logical value indicating whether to compute standard error for the estimates, defaults to FALSE.
#' @param discrepancy a logical value indicating whether to compute discrepancy on the estimates resulting from using saddlepoint likelihood, defaults to FALSE.
#' @param user.ineq.constraint.function a function, user-defined inequality constraint function, defaults to NULL.
#' @param opts.user A named list of options for the nloptr optimizer. This
#'   should be a subset of the following elements: ftol_abs (> 0),
#'   maxeval (a positive integer), xtol_rel (> 0), print_level (0, 1, 2, or 3).
#'   See the nloptr package documentation for more details about these options.
#'   By default, ftol_abs = 0, maxeval = 1e4, xtol_rel = 1.0e-7, print_level = 0.
#'   Note: The algorithm option is fixed to "NLOPT_LD_SLSQP" and cannot be changed by the user.
#' @param sadd.eqn.opts A named list of options for the nloptr optimizer. This sets options for optimization which are
#'   only needed when `std.error = TRUE` or `discrepancy = TRUE`. By default, ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0.
#'   Note: `xtol_rel` argument is set to 1.0e-12 and changing it may affect std.error computations.
#' @return A list of MLEs. If `std.error` is TRUE, the list also includes standard errors of MLEs and inverse Hessian. See TO DO list
#  TO DO : Design this ...
#' @examples
#' \dontrun{
#' # TO DO: Add examples
#' }
#' @export
find.saddlepoint.MLE <- function(observed.data,
                                 model.cgf,
                                 starting.theta,
                                 lb.theta = rep(-Inf, times = length(starting.theta)),
                                 ub.theta = rep(Inf, times = length(starting.theta)),
                                 starting.tvec = rep(0, times = length(observed.data)),
                                 lb.tvec = rep(-Inf, times=length(observed.data)),
                                 ub.tvec = rep(Inf, times=length(observed.data)),
                                 std.error = FALSE,
                                 discrepancy = FALSE,
                                 user.ineq.constraint.function = NULL,
                                 opts.user = list(ftol_abs = 0, maxeval = 1e4, xtol_rel = 1.0e-7, print_level = 0),
                                 sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)
) {

  if (!is(model.cgf, "CGF")) stop("model.cgf must be of class 'CGF'")
  if(length(lb.theta) != length(starting.theta) || length(ub.theta) != length(starting.theta) || !is.numeric(lb.theta) || !is.numeric(ub.theta)) stop("lb.theta or ub.theta has an incorrect length or is not numeric")
  if(length(lb.tvec) != length(starting.tvec) || length(ub.tvec) != length(starting.tvec) || !is.numeric(lb.tvec) || !is.numeric(ub.tvec)) stop("lb.tvec or ub.tvec has an incorrect length or is not numeric")
  if(any(starting.theta < lb.theta) || any(starting.theta > ub.theta)) stop("starting.theta is not within the bounds specified")
  if(!is.numeric(observed.data)) stop("observed.data not defined for ", class(observed.data))
  if (length(starting.tvec) != length(observed.data)) stop("Size of observed.data and starting.tvec arguments do not match")
  
  
  # Combine tvec and theta arguments to a single vector
  a <- c(starting.tvec, starting.theta)
  
  # Define objective, equality constraint and inequality constraint functions
  objective.function <- get.saddlepoint.nll.function(tvec = starting.tvec, theta = starting.theta, model.cgf = model.cgf)
  eq.constraint.function <- get.saddlepoint.eq.constraint.function(tvec = starting.tvec, theta = starting.theta, observed.data = observed.data, model.cgf = model.cgf)
  ineq.constraint.function <- get.ineq.constraint.function(tvec = starting.tvec, theta = starting.theta, model.cgf = model.cgf, user.ineq.constraint.function = user.ineq.constraint.function)
  
  # configure optimizer options
  opts = configure.opts(opts.user) # checks and modifies user-provided options for the optimizer.
  
  # Find the maximum likelihood estimates of the parameters using NLOPT
  MLEs = nloptr::nloptr(x0 = a,
                        eval_f = objective.function,
                        eval_g_eq = eq.constraint.function,
                        eval_g_ineq = ineq.constraint.function,
                        opts = opts,
                        lb =  c(lb.tvec, lb.theta), ub = c(ub.tvec, ub.theta))
  MLEs.tvec = head(MLEs$solution, length(lb.tvec))
  MLEs.theta = tail(MLEs$solution, length(lb.theta))
  
  # Calculate standard errors of MLEs
  if(std.error == TRUE || discrepancy == TRUE){
    res <- compute.std.error(observed.data = observed.data, combined.estimates = MLEs$solution,
                             model.cgf = model.cgf,
                             objective.function = objective.function,
                             # eq.constraint.function = eq.constraint.function,
                             lb.tvec = lb.tvec, ub.tvec = ub.tvec,
                             # ineq.constraint.function = ineq.constraint.function,
                             sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0))
    MLEs$std.error <- res$std.error
    MLEs$inverse.hessian = res$inverse.hessian
  }
  if(discrepancy == TRUE){
    
    grad_FuncT = computeFuncTGradient(tvec = MLEs.tvec,
                                      theta = MLEs.theta,
                                      observations = observed.data,
                                      modelCGF = model.cgf)
    # The formula for the discrepancy should be (-Hessian*gradientOfFuncT) if Hessian is negative definite
    # However, we are minimising the negative log-likelihood which has positive definite Hessian, we therefore
    # omit the negative sign from the formula for this computation
    disc = MLEs$inverse.hessian %*% grad_FuncT$grad.theta.funcT
    MLEs$discrepancy = disc
  }
  
  MLEs$MLEs.tvec = MLEs.tvec
  MLEs$MLEs.theta = MLEs.theta
  MLEs
}


#' @title Find maximum likelihood estimates using the zeroth order saddlepoint approximation.
#' 
#' @description
#' This function uses the nloptr package for optimization to find the maximum likelihood estimates (MLEs) using the zeroth order saddlepoint approximation. It takes the observed data, a CGF object that represents the CGF of the observed data, and various other parameters related to the optimization procedure.
#' @param observed.data a numeric vector - See TO DO list
#' @param model.cgf a CGF object corresponding to the observed data.
#' @param starting.theta a numeric vector of starting values for model parameters.
#' @param lb.theta a numeric vector, vector of lower bounds for model parameters, defaults to -Inf for all parameters.
#' @param ub.theta a numeric vector, vector of upper bounds for model parameters, defaults to Inf for all parameters.
#' @param starting.tvec a numeric vector, vector of starting values for saddlepoint parameters, defaults to zeros.
#' @param lb.tvec a numeric vector, vector of lower bounds for saddlepoint parameters, defaults to -Inf for all parameters.
#' @param ub.tvec a numeric vector, vector of upper bounds for saddlepoint parameters, defaults to Inf for all parameters.
#' @param std.error a logical value indicating whether to compute standard error for the estimates, defaults to FALSE.
#' @param discrepancy a logical value indicating whether to compute discrepancy on the estimates resulting from using saddlepoint likelihood, defaults to FALSE.
#' @param user.ineq.constraint.function a function, user-defined inequality constraint function, defaults to NULL.
#' @param opts.user A named list of options for the nloptr optimizer. 
#' @param sadd.eqn.opts A named list of options for the nloptr optimizer. This sets options for optimization which are
#'  only needed when `std.error = TRUE` or `discrepancy = TRUE`. By default, ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0.
#'  Note: `xtol_rel` argument is set to 1.0e-12 and changing it may affect std.error computations.
#' @return A list of MLEs. If `std.error` is TRUE, the list also includes standard errors of MLEs and inverse Hessian. See TO DO list
#  TO DO : Design this ...
#' @examples
#' \dontrun{
#' # TO DO: Add examples
#' }
#' @export
find.zeroth.saddlepoint.MLE <- function(observed.data,
                                        model.cgf,
                                        starting.theta,
                                        lb.theta = rep(-Inf, times = length(starting.theta)),
                                        ub.theta = rep(Inf, times = length(starting.theta)),
                                        starting.tvec = rep(0, times = length(observed.data)),
                                        lb.tvec = rep(-Inf, times=length(observed.data)),
                                        ub.tvec = rep(Inf, times=length(observed.data)),
                                        std.error = FALSE,
                                        discrepancy = FALSE,
                                        user.ineq.constraint.function = NULL,
                                        opts.user = list(ftol_abs = 0, maxeval = 1e4, xtol_rel = 1.0e-7, print_level = 0),
                                        sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)
){
  
  if (!is(model.cgf, "CGF")) stop("model.cgf must be of class 'CGF'")
  if(length(lb.theta) != length(starting.theta) || length(ub.theta) != length(starting.theta) || !is.numeric(lb.theta) || !is.numeric(ub.theta)) stop("lb.theta or ub.theta has an incorrect length or is not numeric")
  if(any(starting.theta < lb.theta) || any(starting.theta > ub.theta)) stop("starting.theta is not within the bounds specified")
  if(!is.numeric(observed.data)) stop("observed.data not defined for ", class(observed.data))
  if(length(starting.tvec) != length(observed.data)) stop("Size of observed.data and starting.tvec arguments do not match")
  
  # Combine tvec and theta arguments to a single vector
  a <- c(starting.tvec, starting.theta)
  
  # Define objective, equality, and inequality constraint functions
  zeroth.saddlepoint.ll.function <- get.zeroth.saddlepoint.ll.function(tvec = starting.tvec, theta = starting.theta, model.cgf = model.cgf)
        objective.function <- function(a){
          obj_and_grad_list = zeroth.saddlepoint.ll.function(a)
          lapply(obj_and_grad_list, function(x) -x)
        }
  eq.constraint.function <- get.saddlepoint.eq.constraint.function(tvec = starting.tvec, theta = starting.theta, observed.data = observed.data, model.cgf = model.cgf)
  ineq.constraint.function <- get.ineq.constraint.function(tvec = starting.tvec, theta = starting.theta, model.cgf = model.cgf, user.ineq.constraint.function = user.ineq.constraint.function)
  
  # configure optimizer options
  opts = configure.opts(opts.user) # checks and modifies user-provided options for the optimizer.
  
  # Find the maximum likelihood estimates of the parameters using NLOPT
  MLEs = nloptr::nloptr(x0 = a,
                        eval_f = objective.function,
                        eval_g_eq = eq.constraint.function,
                        eval_g_ineq = ineq.constraint.function,
                        opts = opts,
                        lb =  c(lb.tvec, lb.theta), ub = c(ub.tvec, ub.theta))
  MLEs.tvec = head(MLEs$solution, length(lb.tvec))
  MLEs.theta = tail(MLEs$solution, length(lb.theta))
  
  #  # Calculate standard errors of MLEs
  #  if(std.error == TRUE || discrepancy == TRUE){
  #    res <- compute.std.error(observed.data = observed.data, combined.estimates = MLEs$solution,
  #                             model.cgf = model.cgf,
  #                             objective.function = objective.function,
  #                             # eq.constraint.function = eq.constraint.function,
  #                             lb.tvec = lb.tvec, ub.tvec = ub.tvec,
  #                             # ineq.constraint.function = ineq.constraint.function,
  #                             sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0))
  #    MLEs$std.error <- res$std.error
  #    MLEs$inverse.hessian = res$inverse.hessian
  # }
  #  if(discrepancy == TRUE){
  #    grad_FuncT = computeFuncTGradient(tvec = NULL,
  #                                      theta = MLEs.theta,
  #                                      observations = observed.data,
  #                                      modelCGF = model.cgf)
  #    # The formula for the discrepancy should be (-Hessian*gradientOfFuncT) if Hessian is negative definite
  #    # However, we are minimising the negative log-likelihood which has positive definite Hessian, we therefore
  #    # omit the negative sign from the formula for this computation
  #    disc = res$inverse.hessian %*% grad_FuncT$grad.theta.funcT
  #    MLEs$discrepancy = disc
  #  }
  
  MLEs$MLEs.tvec = MLEs.tvec
  MLEs$MLEs.theta = MLEs.theta
  MLEs
}



#' @noRd
configure.opts <- function(opts.user) {
  # default optimizer options
  opts.default = list(algorithm = "NLOPT_LD_SLSQP",
                      ftol_abs = 0,
                      maxeval = 1e3,
                      xtol_rel = 1.0e-12,
                      print_level = 0)
  
  valid.option.names = setdiff(names(opts.default), "algorithm") # Valid option names (excluding "algorithm")
  if (any(names(opts.user) == "")) stop("All elements in opts.user must have names. Valid options are: ", paste(valid.option.names, collapse = ", "))
  if (!all(names(opts.user) %in% valid.option.names)) stop("Invalid option name(s) provided. Valid options are: ", paste(valid.option.names, collapse = ", "))
  modifyList(opts.default, opts.user)
}
