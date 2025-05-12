# R/find-saddlepoint-mle.R
# Objects: find.saddlepoint.MLE


#' @title Find maximum likelihood estimates using the saddlepoint likelihood.
#'
#' @description
#' This function uses the [`nloptr`] package for optimization to find the maximum likelihood estimates (MLEs). 
#' It takes a CGF object along with various parameters related to the optimization procedure. 
#' An alternative zeroth-order method can be enabled by setting the `zeroth.order` parameter to `TRUE`. 
#' 
#' 
#' @details
#' - **Observed Data**: `observed.data` can be a numeric vector or a matrix/data frame.
#'   If it's a matrix/data frame, columns are treated as i.i.d. replicate blocks.
#'   
#' - **Zeroth-order**: The `zeroth.order` argument allows toggling between the standard and zeroth-order saddlepoint likelihoods. 
#'   When `zeroth.order` is set to TRUE, the objective function is modified to the zeroth-order saddlepoint approximation to the negative log-likelihood 
#'   \deqn{\hat{t} y - K_Y (\hat{t}; \theta).} 
#'   The estimates of \eqn{\theta} are calculated based on observations \eqn{y}, 
#'   where \eqn{K_Y(t; \theta)} is the CGF of the underlying random variable, and \eqn{\hat{t}} is the solution 
#'   to the saddlepoint equation. When FALSE, the standard saddlepoint likelihood is used.
#' 
#' 
#' @param observed.data A numeric vector, matrix, or data frame containing the observed data.
#' @param cgf A CGF object corresponding to the distribution of the observed data.
#' @param starting.theta A numeric vector of starting values for model parameters.
#' @param lb.theta A numeric vector of lower bounds for `theta`. Defaults to \code{-Inf}.
#' @param ub.theta A numeric vector of upper bounds for `theta`. Defaults to \code{Inf}.
#' @param starting.tvec A numeric vector of starting values for the saddlepoint \eqn{t}. Defaults to \code{0}.
#' @param lb.tvec A numeric vector, vector of lower bounds for saddlepoint parameters, defaults to -Inf for all parameters.
#' @param ub.tvec A numeric vector, vector of upper bounds for saddlepoint parameters, defaults to Inf for all parameters.
#' @param std.error Logical. If `TRUE`, standard errors of the MLEs are computed. Defaults to `FALSE`.
#' @param discrepancy Logical. If `TRUE`, compute a discrepancy measure that approximates 
#'   how different the saddlepoint approximation-based estimates are from the unknown true MLEs. Defaults to `FALSE`.
#' @param user.ineq.constraint.function Optional user-defined inequality constraints. Defaults to `NULL`.
#' @param opts.user A named list of options for the [`nloptr`] optimizer. This
#'   should be a subset of the following elements: ftol_abs (> 0),
#'   maxeval (a positive integer), xtol_rel (> 0), print_level (0, 1, 2, or 3).
#'   See the [`nloptr`] package documentation for more details about these options.
#'   By default, ftol_abs = 0, maxeval = 1e4, xtol_rel = 1.0e-7, print_level = 0.
#'   Note: The algorithm option is fixed to "NLOPT_LD_SLSQP" and cannot be changed by the user.
#' @param zeroth.order A logical value indicating whether to use the zeroth-order saddlepoint likelihood. Defaults to FALSE, using the standard saddlepoint likelihood.
#' 
#' 
#' @return A list containing:
#'   \itemize{
#'     \item \code{MLEs.tvec, MLEs.theta}: The estimated saddlepoint \eqn{\hat{t}} and parameter vector \eqn{\hat{\theta}}.
#'     \item \code{std.error}, \code{inverse.hessian}: If \code{std.error=TRUE}, these hold the standard errors and the inverse Hessian.
#'     \item \code{discrepancy}: If \code{discrepancy=TRUE}, this is the approximated discrepancy between the resulting saddlepoint approximation-based MLEs and the unknown true MLEs.
#'     \item \code{solution, status, message}: The raw output from [`nloptr`].
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1); x = rgamma(50, shape = 10, rate = 0.5)
#' res <- find.saddlepoint.MLE(observed.data = x, cgf = GammaCGF, starting.theta = c(1,1), std.error=TRUE)
#' res$MLEs.theta
#' res$std.error
#' }
#' @export
find.saddlepoint.MLE <- function(observed.data,
                                 cgf,
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
                                 zeroth.order = FALSE  # Logical flag to use zeroth order method
) {
  
  if (is.matrix(observed.data) || is.data.frame(observed.data)) {
    message("Treating columns of 'observed.data' as i.i.d. replicate blocks.")
    observed.data <- as.numeric(observed.data)
  }
  if (!is(cgf, "CGF")) stop("cgf must be of class 'CGF'")
  if(length(lb.theta) != length(starting.theta) || length(ub.theta) != length(starting.theta) || !is.numeric(lb.theta) || !is.numeric(ub.theta)) stop("lb.theta or ub.theta has an incorrect length or is not numeric")
  if(length(lb.tvec) != length(starting.tvec) || length(ub.tvec) != length(starting.tvec) || !is.numeric(lb.tvec) || !is.numeric(ub.tvec)) stop("lb.tvec or ub.tvec has an incorrect length or is not numeric")
  if(any(starting.theta < lb.theta) || any(starting.theta > ub.theta)) stop("starting.theta is not within the bounds specified")
  if(!is.numeric(observed.data)) stop("observed.data not defined for type ", class(observed.data))
  if (length(starting.tvec) != length(observed.data)) stop("Size of observed.data and starting.tvec arguments do not match")
  
  # Combine tvec and theta arguments into a single vector for the optimizer
  initial_x <- c(starting.tvec, starting.theta)
  
  # Define objective, equality constraint and inequality constraint functions
  objective.function <- get.saddlepoint.nll.function(tvec = starting.tvec, theta = starting.theta, cgf = cgf)
  eq.constraint.function <- get.saddlepoint.eq.constraint.function(tvec = starting.tvec, theta = starting.theta, observed.data = observed.data, cgf = cgf)
  ineq.constraint.function <- get.ineq.constraint.function(tvec = starting.tvec, theta = starting.theta, cgf = cgf, user.ineq.constraint.function = user.ineq.constraint.function)
  
  # Override objective function if zeroth order method is selected
  if(zeroth.order) objective.function <- get.zeroth.saddlepoint.nll.function(tvec = starting.tvec, 
                                                                             theta = starting.theta, 
                                                                             cgf = cgf)
  
  
  # configure optimizer options
  opts = configure.opts(opts.user) # checks and modifies user-provided options for the optimizer.
  
  # Find the maximum likelihood estimates of the parameters using NLOPT
  MLEs = nloptr::nloptr(x0          = initial_x,
                        eval_f      = objective.function,
                        eval_g_eq   = eq.constraint.function,
                        eval_g_ineq = ineq.constraint.function,
                        opts        = opts,
                        lb          = c(lb.tvec, lb.theta), 
                        ub          = c(ub.tvec, ub.theta))
  
  
  MLEs.tvec = head(MLEs$solution, length(lb.tvec))
  MLEs.theta = tail(MLEs$solution, length(lb.theta))
  
  if(MLEs$status < 0){
    cat("Estimates: ", MLEs.theta, "\n")
    stop("Optimization failed: ", MLEs$message)
  } 
  
  # Calculate standard errors of MLEs
  if(std.error || discrepancy){
    res <- compute.std.error(observed.data    = observed.data, 
                             estimated.tvec   = MLEs.tvec,
                             estimated.theta  = MLEs.theta,
                             cgf              = cgf,
                             zeroth.order     = zeroth.order,
                             non.saddlepoint.negll.function = NULL
    )
    MLEs$std.error <- res$std.error
    MLEs$inverse.hessian = res$inverse.hessian
  }
  
  if (discrepancy) {
    # We'll use 'compute.saddlepointLL.correction' with gradient=TRUE to get the derivative 
    # wrt 'theta' of the correction term at the final MLEs
    
    if(zeroth.order) {spa_method <- "zeroth"} else spa_method <- "standard"
    out_ <- compute.saddlepointLL.correction(
        parameter_vector = MLEs.theta,
        observed.data    = observed.data,
        cgf              = cgf,
        tvec.hat = MLEs.tvec,
        gradient = TRUE,     # needed for discrepancy
        hessian  = FALSE,    # not needed
        spa_method   = spa_method
    )
    
    # out_$gradient is the gradient wrt 'theta'
    # The formula for the discrepancy approximation should be (-Hessian*gradientOfFuncT) if Hessian is negative definite
    # However, we are minimising the negative log-likelihood which has positive definite Hessian, we therefore
    # omit the negative sign from the formula for this computation
    discr <- MLEs$inverse.hessian %*% out_$gradient
    MLEs$discrepancy <- as.vector(discr)   
  }
  
  MLEs$MLEs.tvec  = MLEs.tvec
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
