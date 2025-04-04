# R/optimization-functions.R
# Objects: get.saddlepoint.nll.function





#' @title Create the saddlepoint negative log-likelihood function
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
get.saddlepoint.nll.function <- function(tvec, theta, cgf
                                         # , observations) 
                                         ){
  stopifnot(is.numeric(tvec), is.numeric(theta), is(cgf, "CGF"))
  a_init <- c(tvec, theta)
  
  neg_ll <- cgf$.get_private_method("neg_ll")
  neg_ll_spa_wrapper <- function(a){
    tvec_extracted <- a[1:length(tvec)]
    theta_extracted <- a[(length(tvec)+1):length(a)]
    # tvec.hat <- tvec_hat(tvec = tvec_extracted, theta = theta_extracted, 
    #                      observations = observations, K1_fn = cgf$K1, K2_fn = cgf$K2)
    neg_ll(tvec = tvec_extracted, parameter_vector = theta_extracted)
  }
  
  # Build the Tape
  spa_neg_ll_tape <- MakeTape(f = neg_ll_spa_wrapper, x = a_init)
  
  # Return a function that, given a new 'a', returns objective & gradient
  saddlepoint.nll <- function(a) {
    obj_ <- spa_neg_ll_tape(a)  
    gr  <- spa_neg_ll_tape$jacobian(a) 
    list(objective = obj_, gradient = as.vector(gr))
  }
  
}




#' @title Create the zeroth-order saddlepoint negative log-likelihood function
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
get.zeroth.saddlepoint.nll.function <- function(tvec, theta, cgf) {
  stopifnot(is.numeric(tvec), is.numeric(theta), is(cgf, "CGF"))
  a_init <- c(tvec, theta)

  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  neg_tilting_exponent_wrapper <- function(a){
    -tilting_exponent(tvec = a[1:length(tvec)],
                       parameter_vector = tail(a, length(theta)) )
  }
  zeroth_spa_neg_ll_tape <- MakeTape(f = neg_tilting_exponent_wrapper, x = a_init)

  zeroth.saddlepoint.nll <- function(a) {
    obj_ <- zeroth_spa_neg_ll_tape(a)
    gr  <- zeroth_spa_neg_ll_tape$jacobian(a)
    list(objective = obj_, gradient = as.vector(gr))
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
#' 
#' @importFrom Matrix head
#' @importFrom Matrix tail
#' 
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
  
  K1_fn_wrapper <- function(a) cgf$K1(tvec = head(a, length(tvec)), parameter_vector = tail(a, length(theta)))
  K1_fn_tape <- MakeTape(f = K1_fn_wrapper, x = c(tvec, theta))
  
  saddlepoint.eq.constraint.function <- function(a){
    list(constraints = K1_fn_tape(a) - observed.data,
         jacobian =  K1_fn_tape$jacobian(a) )
  }
}





#' @title Create the inequality constraint function
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
        jacobian.with.tvec <- cbind( matrix(0, nrow = nrow(user.ineq.constraint$jacobian), ncol = length(t.vec)),
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
        jacobian.with.tvec <- cbind( matrix(0, 
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



#' @noRd
create_saddlepoint.ineq.constraint_function <- function(tvec, theta, cgf){
  # This is called when the arguments have been checked and are valid
  ineq_fn_wrapper <- function(a) cgf$ineq_constraint(tvec = head(a, length(tvec)), parameter_vector = tail(a, length(theta)))
  ineq_fn_tape <- MakeTape(f = ineq_fn_wrapper, x = c(tvec, theta))
  saddlepoint.ineq.constraint.function <- function(a) {
    list(constraints = ineq_fn_tape(a), jacobian =  ineq_fn_tape$jacobian(a) )
  }
}




#' @title Numerical saddlepoint equation solver
#' @description
#' Solves the saddlepoint equation \eqn{K1(t; \theta) = y} for a given CGF object.
#' It uses numeric optimization (via \code{nloptr}) to find the vector \eqn{t} 
#' that satisfies \eqn{K1(t; \theta) = y}.
#' 
#' If the CGF object offers an analytical solution 
#' (see \code{cgf_object$analytic_tvec_hat} if available), you may prefer that 
#' approach instead.
#' 
#'
#' @importFrom nloptr nloptr
#'
#' @param theta Numeric vector of model parameters for which the CGF is defined.
#' @param y Numeric vector of observations.
#' @param cgf A CGF object (class `"CGF"`).
#' @param starting.tvec Numeric start values for \eqn{t}, defaults to \code{0} for each element.
#' @param lb,ub Numeric vectors for lower and upper bounds of \eqn{t}, 
#'   each having the same length as \code{starting.tvec}. Defaults are \code{-Inf} and \code{Inf}.
#' @param sadd.eqn.opts List of options for the \code{nloptr} optimizer. 
#'   Defaults to \code{list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0)},
#'   with algorithm fixed at \code{"NLOPT_LD_SLSQP"}.
#' @param warn_residual Logical. If \code{TRUE}, compute the residual 
#'   \eqn{\|\mathit{K1}(t^*;\theta)-y\|_\infty} and warn if it exceeds \code{tol}.
#' @param tol Numeric tolerance for residual checks (default \code{1e-4}).
#'   
#' @details
#' Minimizes the function \eqn{K(t; \theta) - \sum(t_i y_i)} to enforce \eqn{K1(t, \theta) = y}.
#'
#' @return A numeric vector \eqn{t} solving \eqn{K1(t, \theta) = y}.
#' 
#' @references
#' \itemize{
#'   \item \href{https://nlopt.readthedocs.io/en/latest/}{NLOpt Documentation}
#' }
#'
#' @examples
#' \dontrun{
#' # TO DO: Add examples
#' }
#'
#' @export
saddlepoint.solve <- function(theta, y, cgf,
                              starting.tvec = rep(0, times = length(y)),
                              lb = rep(-Inf, times=length(y)),
                              ub = rep(Inf, times=length(y)),
                              sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0),
                              warn_residual = TRUE, 
                              tol = 1e-4
                              ) {
  
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  
  if ( cgf$has_analytic_tvec_hat()  ) {
    message("An analytical solution is available via `cgf$analytic_tvec_hat(...)`. ",
            "Consider using that instead of numeric optimization.")
  }
  
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
  
  res = nloptr::nloptr(x0 = starting.tvec,
                       eval_f = eval_f,
                       eval_g_ineq = ineq.as.a.function.of.tvec,
                       opts = opts,
                       lb = lb, ub = ub)
  if (res$status < 0) warning("Saddlepoint solver failed to converge. Consider checking initial values or constraints.")
  
  t_star <- res$solution
  
  if (warn_residual) {
    residual_value <- max(abs(cgf$K1(t_star, theta) - y))
    if (residual_value > tol) {
      warning(
        sprintf("Saddlepoint solution residual = %g exceeds tolerance of %g. ",
                residual_value, tol),
        "Optimization may not have converged to a sufficiently accurate point."
      )
    }
  }
  
  t_star
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
#' @importFrom methods is
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
  method_ <- ifelse(zeroth.order, "zeroth", "standard")
  
  matrix.H <- compute.spa.negll(parameter_vector = estimated.theta, 
                                observed.data    = observed.data,
                                cgf              = cgf, 
                                tvec.hat         = estimated.tvec, 
                                hessian          = TRUE,
                                spa_method       = method_)$hessian
  
  if (!is.null(non.saddlepoint.negll.function)) {
    if (!is.function(non.saddlepoint.negll.function)) stop("'non.saddlepoint.negll.function' must be a function")
    if (!all(c("objective", "hessian") %in% names(non.saddlepoint.negll.function(estimated.theta)))) stop("'non.saddlepoint.negll.function' must have 'objective' and 'hessian' as output")
    matrix.H <- matrix.H + non.saddlepoint.negll.function(estimated.theta)$hessian
  }
  
  inverse.hessian <- solve(matrix.H)
  list(std.error = sqrt(diag(inverse.hessian)), inverse.hessian = inverse.hessian)
}












