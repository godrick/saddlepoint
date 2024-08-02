#' @title Create a saddlepoint negative log-likelihood function
#' @description This function creates and returns a function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param cgf An object of class 'CGF'.
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
get.saddlepoint.nll.function <- function(tvec, theta, cgf){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(cgf, "CGF"))
  adf.negll = makeADFunNegll(tvec = tvec, theta = theta, cgf = cgf$get_ptr())
  
  saddlepoint.nll = function(a){computeCombinedGradient(combined_vector = a, adf = adf.negll)}
}


#' @title Create a zeroth-order saddlepoint negative log-likelihood function
#' @description This function creates and returns a function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param cgf An object of class 'CGF'.
#' @return A function that takes a vector 'a' as an argument. When `a = c(tvec, theta)` is passed to the returned function, it yields a list in the format: \code{list(objective = , gradient = )}.
#  the gradient of the function with respect to both \code{tvec} and \code{theta}.
#' 
#' @examples
#' \dontrun{
#'   TO DO: write a working example
#' }
#' 
#' @export
get.zeroth.saddlepoint.nll.function <- function(tvec, theta, cgf){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(cgf, "CGF"))
  adf.zeroth.nll = makeADFunZerothNegll(tvec = tvec, theta = theta, cgf = cgf$get_ptr())
  
  zeroth.saddlepoint.nll = function(a){
    computeCombinedGradient(combined_vector = a, adf = adf.zeroth.nll)
  }
}









#' @title Create the saddlepoint equality constraint function
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
#' @param cgf An object of class 'CGF'.
#' @return A function that accepts a vector 'a' as an argument. When `a = c(tvec, theta)` is passed to this function, it generates a list containing 'constraints' and 'jacobian'. 'constraints' are computed as \eqn{K'(t;\theta) - y}, and 'jacobian' represents the gradient of these constraints with respect to both \code{tvec} and \code{theta}.
#' @examples
#' \dontrun{
#' TO DO: write a working example
#'   f <- get.saddlepoint_eq_constraint.function(tvec, theta, observed.data, cgf)
#'   f(c(tvec, theta)) # returns a list of the form list(constraints = , jacobian = )
#' }
#' @export
get.saddlepoint.eq.constraint.function <- function(tvec, theta, observed.data, cgf){
  stopifnot(is.numeric(tvec), is.numeric(theta), is.numeric(observed.data), is(cgf, "CGF"))
  adf.K1 = makeADFunK1(tvec = tvec, theta = theta, cgf = cgf$get_ptr())
  
  saddlepoint.eq.constraint.function <- function(a){
    K1.and.grad = computeCombinedGradient(combined_vector = a, adf = adf.K1)
    list(constraints = K1.and.grad$objective - observed.data,
         jacobian = matrix(K1.and.grad$gradient, nrow = length(tvec), byrow = TRUE))
  }
}

#' @title Create an inequality constraint function
#' @description This function wraps and returns an inequality constraint function of the form \code{function(a) \{...\}}, where 'a' combines \code{tvec} and \code{theta} arguments to a single vector, in that order.
#'
#'
#' @details
#' This function constructs and integrates an inequality constraint based on \code{cgf}
#' with any user-defined inequality constraint function specified in \code{user.ineq.constraint.function}.
#' If \code{user.ineq.constraint.function} is `NULL`, the function checks whether the \code{cgf} itself
#' includes an inherent inequality constraint.
#'
#'
#' @param tvec A numeric vector.
#' @param theta A numeric vector.
#' @param cgf An object of class 'CGF'.
#' @param user.ineq.constraint.function An optional additional inequality added by a user. Default is NULL. See TO DO list.
#  TO DO: Add additional documentation on this.
#' @return A function that takes a vector `a` as an argument. This function returns either NULL or a list with 'constraints' and 'jacobian'.
#' @examples
#' \dontrun{
#' TO DO: Add a working example
#'   f <- get.ineq_constraint.function(tvec, theta, cgf)
#'   f(c(tvec, theta)) # returns a list of the form list(constraints = , jacobian = ) or NULL
#' }
#' @export
get.ineq.constraint.function <- function(tvec, theta, cgf, user.ineq.constraint.function = NULL){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(cgf, "CGF"))
  #...
  # This function returns ineq.constraint.function which is either NULL or a function with a single argument 'a'
  # As a result, the output of this function can be directly used in the optimiser.
  m = length(tvec)
  
  
  # First check if the cgf is subject to any constraints, i.e., tvec is constrained
  # We do this by checking if the saddlepoint-based ineq_constraint returns any value
  vector.of.ineq.constraint.values <- cgf$ineq_constraint(tvec, theta)
  
  
  if (!length(vector.of.ineq.constraint.values)) {
    # tvec is not constrained
    # Therefore, for optimisation, ineq.constraint.function will be null
    ineq.constraint.function <- tvec.ineq.constraint.function <- NULL
    
    # If a user defines some constraint on theta, we incorporate the tvec-related jacobian, which is expected by the optimiser
    if (!is.null(user.ineq.constraint.function)) {
      # some checks
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
                                                                                        cgf = cgf)
    
    # If user does not define any constraint on theta, the optimiser only needs the saddlepoint-based constraints (cgf/tvec constraints)
    ineq.constraint.function <- function(a) {
      # Evaluate and return saddlepoint-based inequality constraint function for tvec and theta
      saddlepoint.ineq.constraint.function(a)
    }
    tvec.ineq.constraint.function <- ineq.constraint.function
    
    if (!is.null(user.ineq.constraint.function)) {
      # some checks
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
        jacobian.with.tvec <- cbind(matrix(0, 
                                           nrow = nrow(user.ineq.constraint$jacobian), 
                                           ncol = length(t.vec)),
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
  
  if (is.null(ineq.constraint.function)) { return(NULL) }
  attributes(ineq.constraint.function) <- list(tvec.ineq.constraint.function = tvec.ineq.constraint.function)
  ineq.constraint.function
  
  # attributes(ineq.constraint.function) <- list(tvec.ineq.constraint.function = tvec.ineq.constraint.function)
  # if(!length(ineq.constraint.function)) ineq.constraint.function = NULL
  # ineq.constraint.function
}

#' Create saddlepoint (CGF-based) inequality constraint function
#'
#' This function constructs an inequality constraint function based on the model CGF.
#'
#' @param tvec A numeric vector.
#' @param theta A numeric vector of model parameters.
#' @param cgf An object of class 'CGF'.
#'
#' @return A function that accepts a single vector argument 'a'. When `a = c(tvec, theta)` is passed, the function yields a list in the form \code{list(constraints = , jacobian = )}, where 'constraints' are the calculated inequalities and 'jacobian' is the gradient of the constraints with respect to both \code{tvec} and \code{theta}.
#'
#' @noRd
create_saddlepoint.ineq.constraint_function <- function(tvec, theta, cgf){
  adf.ineq = makeADFunIneqConstraint(tvec = tvec, theta = theta, cgf = cgf$get_ptr())
  
  saddlepoint.ineq.constraint.function <- function(a) {
    ineqConst = computeCombinedGradient(combined_vector = a, adf = adf.ineq)
    list(constraints = ineqConst$objective,
         jacobian = matrix(ineqConst$gradient, nrow = length(ineqConst$objective), byrow = TRUE))
  }
}






#' Solution of the saddlepoint equation
#'
#' This function minimises the expression \eqn{K'(t, \theta)-y} of the saddlepoint equation, and returns
#' the value of \eqn{t} for which \eqn{K'(t, \theta) = y}, for a specified/known value of \eqn{\theta}.
#'
#' @importFrom nloptr nloptr
#'
#' @param theta A numeric vector of model parameters for which the CGF is defined.
#' @param y A numeric vector of observations.
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
saddlepoint.eqn.solve <- function(theta, y, cgf,
                                  starting.tvec = rep(0, times = length(y)),
                                  lb = rep(-Inf, times=length(y)),
                                  ub = rep(Inf, times=length(y)),
                                  sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)) {
  
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  if(length(lb) != length(starting.tvec) || length(ub) != length(starting.tvec) || !is.numeric(lb) || !is.numeric(ub)) stop("lb or ub has an incorrect length or is not numeric")

  # get tvec-related ineq.constraint_function
  tvec.ineq.constraint.function = attributes(get.ineq.constraint.function(tvec = starting.tvec, theta = theta,
                                                                          cgf = cgf,
                                                                          user.ineq.constraint.function = NULL))$tvec.ineq.constraint.function
  objective.fun <- function(t.vec){
    cgf$K(t.vec, theta) - sum(t.vec*y)
  }
  grad.objective.fun <- function(t.vec){
    cgf$K1(t.vec, theta) - y
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
#' This function calculates the standard error and the inverse of the Hessian matrix 
#' for a model incorporating both the saddlepoint likelihood and an optional 
#' non-saddlepoint negative log-likelihood component.
#' 
#'
#' @param observed.data A numeric vector of observations used in the saddlepoint likelihood calculations.
#' @param estimated.tvec A numeric vector of MLE values for the saddlepoint values `tvec`.
#' @param estimated.theta A numeric vector of MLE values for the model parameters.
#' @param cgf A CGF object used in the saddlepoint likelihood computation.
#' @param zeroth.order A logical value indicating whether the zeroth-order saddlepoint likelihood is used. Default is FALSE.
#' @param non.saddlepoint.negll.function An optional function specifying a non-saddlepoint negative log-likelihood. If your model used for estimation incorporates a likelihood component that is not based on the saddlepoint approximation, this function must be provided. See details for more information.
#'
#' @details
#' The `non.saddlepoint.negll.function`, if provided, should have the form:
#' \deqn{function(theta) \{ return(list(objective = ..., hessian = ...)) \}}
#' It adds flexibility by allowing incorporation of likelihood components not based on the saddlepoint likelihood.
#' It is important when the analysis combines the saddlepoint likelihood with an external likelihood. 
#' If not provided, the function defaults to using solely the saddlepoint likelihood.
#'
#' @return A list containing:
#'   \itemize{
#'     \item std.error: Standard error computed from the diagonal of the inverse Hessian.
#'     \item inverse.hessian: Inverse of the Hessian matrix.
#'   }
#'
#' @examples
#' # TODO: Add examples
#'
#' @export
compute.std.error <- function(observed.data, 
                              estimated.tvec,
                              estimated.theta,
                              cgf, 
                              zeroth.order = FALSE,
                              non.saddlepoint.negll.function = NULL) {
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  computeNegll_fn <- if (zeroth.order) computeZerothNegll else computeNegll
  matrix.H <- matrix(computeNegll_fn(tvec = estimated.tvec, theta = estimated.theta,
                                    observations = observed.data, cgf = cgf$get_ptr())$hessian,
                    nrow = length(estimated.theta))

  
  if (!is.null(non.saddlepoint.negll.function)) {
    if (!is.function(non.saddlepoint.negll.function)) stop("'non.saddlepoint.negll.function' must be a function")
    matrix.H <- matrix.H + non.saddlepoint.negll.function(estimated.theta)$hessian
  }
  
  inverse.hessian <- solve(matrix.H)
  list(std.error = sqrt(diag(inverse.hessian)),
       inverse.hessian = inverse.hessian)
}



















