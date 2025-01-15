# R/GammaCGF.R
# Objects: GammaCGF, GammaModelCGF

#' Gamma CGF Object
#'
#' A ready-to-use CGF object for the Gamma distribution with shape \eqn{\alpha}
#' and rate \eqn{\beta}. The \code{parameter_vector} used when calling methods like `K(tvec, parameter_vector)`
#' shoulde be a numeric vector \eqn{c(\alpha, \beta)}.
#'
#' @details
#' **CGF**: For a Gamma random variable \eqn{X} with shape \eqn{\alpha} and rate
#' \eqn{\beta}, the cumulant generating function is:
#' \deqn{K(t;\alpha, \beta) = -\alpha \,\log \bigl(1 - t/\beta\bigr), \quad t < \beta.}
#'
#' **Parameter Vector**: The \code{parameter_vector} is assumed to have the form 
#' \eqn{(\alpha, \beta)}. You must ensure 
#' that \eqn{t < \beta} for all evaluations. This object enforces that constraint 
#' as an inequality \eqn{t - \beta < 0}.
#'
#' @format An object of class \code{CGF} (an R6 class), with the usual methods:
#' \code{K, K1, K2, K3operator, K4operator}, etc.
#'
#' @examples
#' # Evaluate K at t=0.5 for shape=2, rate=2 (thus t<2).
#' # param = c(2, 2)
#' GammaCGF$K(0.5, c(2,2))
#'
#' @export
GammaCGF <- createCGF_fromVectorisedFunctions(
  
  K_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K(t) = -alpha * log1p(-t/beta)
    # using log1p(x) for x = -(tvec/beta)
    -alpha * log1p(-tvec/beta)
  },
  
  ########################
  ## K1(t)
  ########################
  K1_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K1(t) = alpha / (beta - t)
    alpha / (beta - tvec)
  },
  
  
  K2_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    # K2(t) = alpha / (beta - t)^2
    alpha / (beta - tvec)^2
  },
  
  
  K3_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    2 * alpha / (beta - tvec)^3
  },
  
  K4_vectorized_func = function(tvec, param) {
    alpha <- param[1]
    beta  <- param[2]
    6 * alpha / (beta - tvec)^4
  },
  
  
  # We want tvec < beta => tvec - beta < 0
  # We'll just return (tvec - beta), which must be negative
  ineq_constraint_func = function(tvec, param) {
    beta <- param[2]
    tvec - beta
  },
  
  analytic_tvec_hat_func = function(y, param)  { param[2] - (param[1] / y) },
  
  param_adaptor = function(x) x[1:2],
  op_name       = "GammaCGF"
)


#' Create a Parametric Gamma CGF Object
#'
#' @description
#' Creates a CGF object for the Gamma distribution with shape \eqn{\alpha(\theta)} and 
#' rate \eqn{\beta(\theta)} defined by user-provided parameter functions. 
#' The user supplies a function that maps their full parameter vector 
#' \code{theta} to \eqn{(\alpha, \beta)}.
#'
#' @param shape A function that accepts a single parameter vector \code{theta} and returns the shape parameter.
#' @param rate A function that accepts a single parameter vector \code{theta} and returns the rate parameter.
#' @param iidReps A positive integer specifying the number of i.i.d. replicates. Defaults to \code{1}.
#' @param ... Additional arguments passed to the underlying CGF creation function. 
#'   This can include optional overrides to the base CGF class.
#'
#'
#' @return A CGF object.
#' @export
GammaModelCGF <- function(shape,
                          rate,
                          iidReps = 1,
                          ...) {
  shape_fn <- validate_function_or_adaptor(shape)
  rate_fn <- validate_function_or_adaptor(rate)
  
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1 || iidReps != as.integer(iidReps)) stop("'iidReps' must be a positive integer.")
  check_tvec <- function(tvec) {
    if (length(tvec) != iidReps) stop(sprintf("`tvec` must have length %d (got %d).", iidReps, length(tvec)))
  }
  
  
  
  # param_adaptor => shape_rate
  base_cgf <- createCGF_fromVectorisedFunctions(
    K_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      -alpha * log1p(-tvec / beta)
    },
    K1_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      alpha / (beta - tvec)
    },
    K2_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      alpha / (beta - tvec)^2
    },
    K3_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      2*alpha / (beta - tvec)^3
    },
    K4_vectorized_func = function(tvec, sr) {
      check_tvec(tvec)
      alpha <- sr[1]
      beta  <- sr[2]
      6*alpha / (beta - tvec)^4
    },
    ineq_constraint_func = function(tvec, sr) {
      check_tvec(tvec)
      beta <- sr[2]
      tvec - beta
    },
    param_adaptor = function(x) c(shape_fn(x), rate_fn(x)),
    analytic_tvec_hat_func = function(y, sr)  {
      if (length(y) != iidReps) stop(sprintf("`y` must have length %d (got %d).", iidReps, length(y)))
      alpha <- sr[1]
      beta  <- sr[2]
      beta - (alpha / y)
    },
    op_name = "GammaModelCGF",
    ...
  )
  
  base_cgf
}











#' A CGF object for a vector of Non-Identical Gamma variables
#'
#' @description
#' Creates a CGF object for a \eqn{d}-dimensional vector of independent, non-identically distributed 
#' Gamma random variables. Each coordinate \eqn{i} is associated with its own shape \eqn{\alpha_i} 
#' and rate \eqn{\beta_i}, determined by two user-provided functions: 
#' \code{shapes} and \code{rates}. 
#' 
#' The \code{shapes} function maps a user-specified parameter vector (e.g., \eqn{\theta}) 
#' to shape parameters \eqn{(\alpha_1, \ldots, \alpha_d)}, while the \code{rates} function 
#' maps the same \eqn{\theta} to rate parameters \eqn{(\beta_1, \ldots, \beta_d)}.
#'
#' While the focus is on modeling a single vector of independent, non-identical Gamma variables, 
#' the design also supports IID replication of these vectors. The \code{iidReps} argument specifies 
#' how many independent copies of the \eqn{d}-dimensional vector are being modeled. Thus, 
#' each input \code{tvec} provided to CGF methods (such as \code{K}, \code{K1}, etc.) must have 
#' length \eqn{d \times \text{iidReps}}.
#' 
#' 
#' @param shapes A function of one argument (\eqn{\theta}) returning a numeric vector 
#'   \eqn{(\alpha_1, \ldots, \alpha_d)} of shape parameters for each Gamma variable.
#' @param rates A function of one argument (\eqn{\theta}) returning a numeric vector 
#'   \eqn{(\beta_1, \ldots, \beta_d)} of rate parameters for each Gamma variable.
#' @param iidReps Integer. The number of IID replicates of the \eqn{d}-dimensional Gamma vector.
#'   Must be a positive integer. Defaults to \code{1}.
#' @param ... Additional arguments passed to the underlying CGF creation function.
#'
#' @return A CGF object
#'
#' @examples
#' \dontrun{
#' # Suppose we have d=2 distinct Gamma variables in a single vector
#' # and want to replicate it iidReps=3 times for a total dimension of 6.
#' # shapes_fn <- function(theta) { c(theta[1], theta[2]) }  # alpha1, alpha2
#' # rates_fn  <- function(theta) { c(theta[3], theta[4]) }  # beta1,  beta2
#' #
#' # gamma_cgf_obj <- GammaNonIdenticalModelCGF(
#' #    shapes = shapes_fn,
#' #    rates  = rates_fn,
#' #    iidReps=3
#' # )
#' #
#' # # Now tvec must be length 6 
#' # tvec  <- rep(0.1, 6)
#' # param <- c(2.0, 1.5, 5.0, 4.0)  # example
#' #
#' # gamma_cgf_obj$K(tvec, param)    # Evaluate K
#' # gamma_cgf_obj$K1(tvec, param)   # Evaluate first derivative
#' # ...
#' }
#'
#' @export
GammaNonIdenticalModelCGF <- function(shapes, rates, iidReps = 1, ...) {
  
  if (!is.numeric(iidReps) || length(iidReps) != 1 || iidReps < 1 || iidReps != as.integer(iidReps)) stop("`iidReps` must be a positive integer.")
  s_fn <- validate_function_or_adaptor(shapes)
  r_fn <- validate_function_or_adaptor(rates)
  
  # We create our param_adaptor to unify shapes(...) and rates(...) into
  # one concatenated vector: c(alpha_1,...,alpha_d, beta_1,...,beta_d).
  param_adaptor <- function(theta) {
    alpha_vec <- s_fn(theta)
    beta_vec  <- r_fn(theta)
    if (length(alpha_vec) != length(beta_vec)) stop("`shapes(theta)` and `rates(theta)` must return vectors of the same length.")
    c(alpha_vec, beta_vec)
  }
  
  createCGF_fromVectorisedFunctions(
    
    # K(tvec) = sum over i of [ - alpha_i * log(1 - t_i / beta_i ) ]
    # We vectorize, so we do it for alpha_rep[i], beta_rep[i], tvec[i].
    K_vectorized_func = function(tvec, param) {
      len_sr <- length(param) # param is the concatenated shapes/rates from param_adaptor
      if (len_sr %% 2 != 0) stop("The concatenated shapes/rates vector must have even length.")
      d <- len_sr / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).", 
                     length(tvec), d * iidReps))
      }
      
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      
      # replicate them for each i.i.d. copy
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      - alpha_rep * log1p(- tvec / beta_rep)
    },
    
    # K1(tvec) = alpha_i / (beta_i - t_i)
    K1_vectorized_func = function(tvec, param) {
      d <- length(param) / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).",
                     length(tvec), d * iidReps))
      }
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      alpha_rep / (beta_rep - tvec)
    },
    
    # K2(tvec) = alpha_i / (beta_i - t_i)^2
    K2_vectorized_func = function(tvec, param) {
      d <- length(param) / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).",
                     length(tvec), d * iidReps))
      }
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      alpha_rep / (beta_rep - tvec)^2
    },
    
    # K3(tvec) = 2 * alpha_i / (beta_i - t_i)^3
    K3_vectorized_func = function(tvec, param) {
      d <- length(param) / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).",
                     length(tvec), d * iidReps))
      }
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      2 * alpha_rep / (beta_rep - tvec)^3
    },
    
    # K4(tvec) = 6 * alpha_i / (beta_i - t_i)^4
    K4_vectorized_func = function(tvec, param) {
      d <- length(param) / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).",
                     length(tvec), d * iidReps))
      }
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      6 * alpha_rep / (beta_rep - tvec)^4
    },
    
    # ineq_constraint_func: We must have t_i < beta_i for each coordinate
    # so we return (t_i - beta_i) which must be < 0.
    ineq_constraint_func = function(tvec, param) {
      d <- length(param) / 2
      if (length(tvec) != d * iidReps) {
        stop(sprintf("`tvec` has length %d; expected %d (d*iidReps).",
                     length(tvec), d * iidReps))
      }
      beta_vec  <- param[(d+1):(2*d)]
      beta_rep  <- rep(beta_vec, times = iidReps)
      
      tvec - beta_rep
    },
    
    param_adaptor = param_adaptor,
    
    # analytic_tvec_hat_func: For each i,
    # t^*_i = beta_i - (alpha_i / y_i).
    analytic_tvec_hat_func = function(y, param) {
      d <- length(param) / 2
      if (length(y) != d * iidReps) {
        stop(sprintf("`y` must have length %d; got %d.", d * iidReps, length(y)))
      }
      alpha_vec <- param[1:d]
      beta_vec  <- param[(d+1):(2*d)]
      alpha_rep <- rep(alpha_vec, times = iidReps)
      beta_rep  <- rep(beta_vec,  times = iidReps)
      
      beta_rep - (alpha_rep / y)
    },
    
    op_name = "GammaNonIdenticalModelCGF",
    ...
  )
}












