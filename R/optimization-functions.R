#' @title Create a saddlepoint negative log-likelihood function
#' @description This function creates and returns a function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param model.cgf An object of class 'CGF'.
#' @return A function that takes a vector 'a' as an argument. When `a = c(tvec, theta)` is passed to the returned function, it yields a list in the format: \code{list(objective = , gradient = )}.
#  the gradient of the function with respect to both \code{tvec} and \code{theta}.
#'
#' @examples
#' \dontrun{
#'   TO DO: write a working example
#'   f <- get.saddlepoint.nll.function(tvec, theta, PoissonCGF)
#'   f(c(tvec, theta)) # returns a list of the form list(objective = , gradient = )
#' }
#' @export
get.saddlepoint.nll.function <- function(tvec, theta, model.cgf){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(model.cgf, "CGF"))
  ADfun.negll = make.ADFunNegll(tvec = tvec, theta = theta, ModelCGF = model.cgf)
  
  saddlepoint.nll = function(a){
    negll_with_gradient(combined_vector = a, ADfun_negll = ADfun.negll)
  }
}

#' @title Create saddlepoint equality constraint function
#'
#' @description
#' This function creates and returns a function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#'
#' @details
#' In this set up, the expression of the saddlepoint equation is defined by \eqn{K'(t;\theta) - y = 0} for \code{observed.data = y}. The returned function is designed to compute both \eqn{K'(t;\theta) - y} and its gradient with respect to both \code{tvec} and \code{theta}.
#'
#'
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param observed.data A numeric vector. See TO DO list.
#  TO DO: Add something on the length ...
#' @param model.cgf An object of class 'CGF'.
#' @return A function that accepts a vector 'a' as an argument. When `a = c(tvec, theta)` is passed to this function, it generates a list containing 'constraints' and 'jacobian'. 'constraints' are computed as \eqn{K'(t;\theta) - y}, and 'jacobian' represents the gradient of these constraints with respect to both \code{tvec} and \code{theta}.
#' @examples
#' \dontrun{
#' TO DO: write a working example
#'   f <- get.saddlepoint_eq_constraint.function(tvec, theta, observed.data, model.cgf)
#'   f(c(tvec, theta)) # returns a list of the form list(constraints = , jacobian = )
#' }
#' @export
get.saddlepoint.eq.constraint.function <- function(tvec, theta, observed.data, model.cgf){
  stopifnot(is.numeric(tvec), is.numeric(theta), is.numeric(observed.data), is(model.cgf, "CGF"))
  ADfun.K1 = make.ADFunK1(tvec = tvec, theta = theta, ModelCGF = model.cgf)
  
  saddlepoint.eq.constraint.function <- function(a){
    K1.and.grad = K1_with_gradient(combined_vector = a, ADfunK1 = ADfun.K1)
    list(constraints = K1.and.grad$fn - observed.data,
         jacobian = matrix(K1.and.grad$gr, nrow = length(tvec), byrow = TRUE))
  }
}

#' @title Create an inequality constraint function
#' @description This function wraps and returns an inequality constraint function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#'
#'
#' @details
#' This function constructs and integrates an inequality constraint based on \code{model.cgf}
#' with any user-defined inequality constraint function specified in \code{user.ineq.constraint.function}.
#' If \code{user.ineq.constraint.function} is `NULL`, the function checks whether the \code{model.cgf} itself
#' includes an inherent inequality constraint.
#'
#'
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param model.cgf An object of class 'CGF'.
#' @param user.ineq.constraint.function An optional additional inequality added by a user. Default is NULL. See TO DO list.
#  TO DO: Add additional documentation on this.
#' @return A function that takes a vector `a` as an argument. This function returns either NULL or a list with 'constraints' and 'jacobian'.
#' @examples
#' \dontrun{
#' TO DO: Add a working example
#'   f <- get.ineq_constraint.function(tvec, theta, model.cgf)
#'   f(c(tvec, theta)) # returns a list of the form list(constraints = , jacobian = ) or NULL
#' }
#' @export
get.ineq.constraint.function <- function(tvec, theta, model.cgf, user.ineq.constraint.function = NULL){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(model.cgf, "CGF"))
  #...
  # This function returns ineq.constraint.function which is either NULL or a function with a single argument 'a'
  # As a result, the output of this function can be directly used in the optimiser.
  m = length(tvec)
  
  
  # First check if the model.cgf is subject to any constraints, i.e., tvec is constrained
  # We do this by checking if the saddlepoint-based ineq_constraint returns any value
  vector.of.ineq.constraint.values <- ineq_constraint(tvec, theta, model.cgf)
  
  if (!length(vector.of.ineq.constraint.values)) {
    # tvec is not constrained
    # Therefore, for optimisation, ineq.constraint.function will be null
    ineq.constraint.function <- tvec.ineq.constraint.function <- NULL
    
    # If a user defines some constraint on theta, we incorporate the tvec-related jacobian, which is expected by the optimiser
    if (!is.null(user.ineq.constraint.function)) {
      # SOme checks
      if (!is.function(user.ineq.constraint.function)) stop("user.ineq.constraint.function must be a function")
      if (!all(c("constraints", "jacobian") %in% names(user.ineq.constraint.function(theta)))) stop("user.ineq.constraint.function must have 'constraints' and 'jacobian' as output")
      
      # Define the new ineq.constraint.function that includes tvec
      ineq.constraint.function <- function(a) {
        # Extract tvec and theta
        t.vec <- a[1:m]
        theta <- a[-(1:m)]
        # Evaluate user-supplied inequality constraint function for theta
        user.ineq.constraint <- user.ineq.constraint.function(theta)
        
        # Add zeros to the jacobian for the tvec components
        jacobian.with.tvec <- cbind(matrix(0, nrow = nrow(user.ineq.constraint$jacobian), ncol = length(t.vec)),
                                    user.ineq.constraint$jacobian)
        
        # Return the constraints and jacobian including tvec
        list(constraints = user.ineq.constraint$constraints,
             jacobian = jacobian.with.tvec)
      }
    }
  } else { # tvec is constrained
    
    saddlepoint.ineq.constraint.function <- create_saddlepoint.ineq.constraint_function(tvec = tvec,
                                                                                        theta = theta,
                                                                                        model.cgf = model.cgf)
    
    # If user does not define any constraint on theta, the optimiser only needs the saddlepoint-based constraints
    ineq.constraint.function <- function(a) {
      # Evaluate are return saddlepoint-based inequality constraint function for tvec and theta
      saddlepoint.ineq.constraint.function(a)
    }
    tvec.ineq.constraint.function <- ineq.constraint.function
    
    if (!is.null(user.ineq.constraint.function)) {
      # Check if user.ineq.constraint.function is a function
      if (!is.function(user.ineq.constraint.function)) stop("user.ineq.constraint.function must be a function")
      if (!all(c("constraints", "jacobian") %in% names(user.ineq.constraint.function(theta)))) stop("user.ineq.constraint.function must have 'constraints' and 'jacobian' as output")
      
      # Define the new ineq.constraint.function that includes tvec
      ineq.constraint.function <- function(a) {
        t.vec <- a[1:m]
        theta <- a[-(1:m)]
        
        # Evaluate user-supplied inequality constraint function for theta
        user.ineq.constraint <- user.ineq.constraint.function(theta)
        
        # Evaluate saddlepoint-based inequality constraint function for tvec and theta
        saddlepoint.ineq.constraint <- saddlepoint.ineq.constraint.function(a)
        
        # Add zeros to the jacobian for the tvec components
        jacobian.with.tvec <- cbind(matrix(0, nrow = nrow(user.ineq.constraint$jacobian), ncol = length(saddlepoint.ineq.constraint$constraints)),
                                    user.ineq.constraint$jacobian)
        
        # Combine the constraints from both functions
        combined.constraints <- c(saddlepoint.ineq.constraint$constraints, user.ineq.constraint$constraints)
        
        # Combine the jacobians from both functions
        combined.jacobian <- rbind(saddlepoint.ineq.constraint$jacobian, jacobian.with.tvec)
        
        # Return the constraints and jacobian including tvec
        list(constraints = combined.constraints,
             jacobian = combined.jacobian)
      }
    }
    
  }
  attributes(ineq.constraint.function) <- list(tvec.ineq.constraint.function = tvec.ineq.constraint.function)
  if(!length(ineq.constraint.function)) ineq.constraint.function = NULL
  return(ineq.constraint.function)
}

#' Create saddlepoint (CGF-based) inequality constraint function
#'
#' This function constructs an inequality constraint function based on the model CGF.
#'
#' @param tvec A numeric vector.
#' @param theta A numeric vector of model parameters.
#' @param model.cgf An object of class 'CGF'.
#'
#' @return A function that accepts a single vector argument 'a'. When `a = c(tvec, theta)` is passed, the function yields a list in the form \code{list(constraints = , jacobian = )}, where 'constraints' are the calculated inequalities and 'jacobian' is the gradient of the constraints with respect to both \code{tvec} and \code{theta}.
#'
#' @noRd
create_saddlepoint.ineq.constraint_function <- function(tvec, theta, model.cgf){
  ADfun.ineq = make.ADFunIneqConstraint(tvec = tvec, theta = theta, ModelCGF = model.cgf)
  
  saddlepoint.ineq.constraint.function <- function(a) {
    ineqConst = ineqConstraint_with_gradient(combined_vector = a, ADfun_ineqConstraint = ADfun.ineq)
    list(constraints = ineqConst$fn,
         jacobian = matrix(ineqConst$gr, nrow = length(ineqConst$fn), byrow = TRUE))
  }
}



#' Solution of the saddlepoint equation
#'
#' This function minimises the expression \eqn{K'(t, \theta)-y} of the saddlepoint equation, and returns
#' the value of \eqn{t} for which \eqn{K'(t, \theta) = y}, for a specified/known value of \eqn{\theta}.
#'
#' @importFrom nloptr nloptr
#'
#' @param y A numeric vector of observations.
#' @param theta A numeric vector of model parameters for which the CGF is defined.
#' @param cgf An object of type "CGF".
#' @param starting.tvec a numeric vector, vector of starting values for saddlepoint parameters, defaults to zeros.
#' @param lb a numeric vector, vector of lower bounds for saddlepoints tvec, defaults to -Inf.
#' @param ub a numeric vector, vector of upper bounds for saddlepoints tvec, defaults to Inf.
#' @param sadd.eqn.opts A list of options for the nlopt optimizer. This should be a named list, where the names are option names and the list elements are the corresponding option values.
#'    The default options are: \code{list(algorithm="NLOPT_LD_SLSQP", ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)}.
#'    Note that the option "algorithm" cannot be changed. See the documentation for the \code{nloptr} package for more details about the options.
#'
#' @return The optimal value of \code{t.vec} for which K1(t.vec, theta) = y, according to the nlopt optimization.
#'
#' @examples
#' # TO DO: Add examples
#'
#' @export
sadd.eqn.fn <- function(theta, y, cgf,
                        starting.tvec = rep(0, times = length(y)),
                        # tvec.ineq.constraint.function = NULL,
                        lb = rep(-Inf, times=length(y)),
                        ub = rep(Inf, times=length(y)),
                        sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)) {
  
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  if(length(lb) != length(starting.tvec) || length(ub) != length(starting.tvec) || !is.numeric(lb) || !is.numeric(ub)) stop("lb.tvec or ub.tvec has an incorrect length or is not numeric")
  # if(!is.null(tvec.ineq.constraint.function) && !is.function(tvec.ineq.constraint.function)) stop("tvec.ineq.constraint.function is not defined for ", class(tvec.ineq.constraint.function))
  
  # get tvec-related ineq.constraint_function
  tvec.ineq.constraint.function = attributes(get.ineq.constraint.function(tvec = starting.tvec, theta = theta,
                                                                          model.cgf = cgf,
                                                                          user.ineq.constraint.function = NULL))$tvec.ineq.constraint.function
  objective.fun <- function(t.vec){
    val = K(t.vec, theta, cgf) - t(as.matrix(t.vec)) %*% y
    # if(is.na(val)) val = Inf
    val
  }
  grad.objective.fun <- function(t.vec){
    temp = as.vector(K1(t.vec, theta, cgf))
    # temp[is.na(temp)]  = -9999
    temp - y
  }
  eval_f <- function(t.vec) {
    list(objective = objective.fun(t.vec),
         gradient = grad.objective.fun(t.vec))
  }
  
  # configure optimizer options
  # The function checks and modifies user-provided options for the optimizer, ensuring they are valid and complete.
  opts = configure.sadd.eqn.opts(sadd.eqn.opts)
  
  ineq.as.a.function.of.tvec = NULL
  if(!is.null(tvec.ineq.constraint.function)){
    # The argument 'a' of tvec.ineq.constraint.function is a combined function of t.vec and theta
    ineq.as.a.function.of.tvec <- function(t.vec) {
      # Combine t.vec and theta into a single vector
      a <- c(t.vec, theta)
      
      temp.jacobian = tvec.ineq.constraint.function(a)$jacobian
      list(constraints = tvec.ineq.constraint.function(a)$constraints,
           jacobian = temp.jacobian[, which(rep(1:(length(a)), length.out = ncol(temp.jacobian)) <= length(t.vec)), drop = FALSE]
      )
    }
  }
  
  nloptr::nloptr(x0 = starting.tvec,
                 eval_f = eval_f,
                 eval_g_ineq = ineq.as.a.function.of.tvec,
                 opts = opts,
                 lb = lb, ub = ub)$solution
}

#' @importFrom utils modifyList
#' @noRd
configure.sadd.eqn.opts <- function(sadd.eqn.opts) {
  # default optimizer options
  sadd.eqn.opts.default = list(algorithm = "NLOPT_LD_SLSQP",
                               ftol_abs = 0,
                               maxeval = 1e3,
                               xtol_rel = 1.0e-12,
                               print_level = 0)
  
  valid.option.names = setdiff(names(sadd.eqn.opts.default), "algorithm") # Valid option names (excluding "algorithm")
  if (any(names(sadd.eqn.opts) == "")) stop("All elements in sadd.eqn.opts must have names. Valid options are: ", paste(valid.option.names, collapse = ", "))
  if (!all(names(sadd.eqn.opts) %in% valid.option.names)) stop("Invalid option name(s) provided. Valid options are: ", paste(valid.option.names, collapse = ", "))
  modifyList(sadd.eqn.opts.default, sadd.eqn.opts)
}




#' Compute standard error and inverse Hessian
#'
#' This function calculates the standard error and the inverse of the Hessian matrix ....
#'
#' @importFrom numDeriv hessian
#' @importFrom utils head tail
#'
#' @param observed.data A numeric vector of observed data
#' @param combined.estimates TO DO A numeric vector of combined estimated parameters
#' @param model.cgf A CGF object
#' @param objective.function A function that defines the objective to be minimized...TO DO
#' @param lb.tvec A numeric vector specifying the lower bounds for tvec (default: -Inf)
#' @param ub.tvec A numeric vector specifying the upper bounds for tvec (default: Inf)
#' @param sadd.eqn.opts A list of options for the nlopt optimizer. This should be a named list, where the names are option names and the list elements are the corresponding option values.
#'    The default options are: \code{list(algorithm="NLOPT_LD_SLSQP", ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)}.
#'    Note that the option "algorithm" cannot be changed. See the documentation for the \code{nloptr} package for more details about the options.
# #' @param ineq.constraint.function An optional function defining inequality constraints (default: NULL)
#'
#' @return A list containing the standard error (std.error) and the inverse Hessian matrix (inverse.hessian)
#'
#' @examples
#' # TODO: Add examples
#'
#' @export
compute.std.error <- function(observed.data, combined.estimates,
                              model.cgf, objective.function,
                              # eq.constraint.function,
                              lb.tvec = rep(-Inf, times=length(observed.data)),
                              ub.tvec = rep(Inf, times=length(observed.data)),
                              sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)
                              # ineq.constraint.function = NULL
) {
  
  estimated.tvec <- head(combined.estimates, length(observed.data))
  estimated.theta <- tail(combined.estimates, length(combined.estimates) - length(estimated.tvec))
  
  
  # Define a function to evaluate the objective function as a function of theta
  nll.as.a.function.of.theta <- function(theta){
    tvec = sadd.eqn.fn(theta = theta, y = observed.data,
                       cgf = model.cgf, starting.tvec = estimated.tvec,
                       lb = lb.tvec, ub = ub.tvec,
                       sadd.eqn.opts = sadd.eqn.opts)
    objective.function(c(tvec, theta))$objective
  }
  
  matrix.H <- numDeriv::hessian(nll.as.a.function.of.theta, estimated.theta)
  
  inverse.hessian <- solve(matrix.H)
  list(std.error = sqrt(diag(inverse.hessian)),
       inverse.hessian = inverse.hessian)
}



















