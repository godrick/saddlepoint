# R/C-BinomialCGF.R






#' Binomial CGF Object
#'
#' A ready-to-use CGF object for the Binomial distribution. This object expects
#' calls of the form \code{BinomialCGF$K(tvec, param)}, where \code{param}
#' is a numeric vector \eqn{(n, p)}:
#' \describe{
#'   \item{n}{A positive integer (the number of trials).}
#'   \item{p}{A probability in \eqn{(0,1)}.}
#' }
#'
#' By default, \code{BinomialCGF} supports vectorized evaluation for i.i.d. replicates of
#' Binomial random variables.
#'
#'
#' @format An object of class \code{CGF} (R6), with standard methods:
#'   \code{K}, \code{K1}, \code{K2}, \code{K3operator}, etc.
#'
#' @examples
#' # Evaluate the expected value of X ~ Binomial(10, 0.3) using the BinomialCGF object
#' BinomialCGF$K1(0, c(10, 0.3))
#'
#' @export
BinomialCGF <- createCGF_fromVectorisedFunctions(
  K_vectorized_func = function(tvec, param) {
    # param = c(n, p)
    # If length(tvec) > 1 => sum_{i=1..length(tvec)} K( tvec[i], param )
    n_val <- param[1]
    p_val <- param[2]
    n_val*log(1 - p_val + p_val*exp(tvec)) 
  },
  K1_vectorized_func = function(tvec, param) {
    n_val <- param[1]
    p_val <- param[2]
    n_val*p_val*exp(tvec) / (1 - p_val + p_val*exp(tvec))
  },
  K2_vectorized_func = function(tvec, param) {
    n_val <- param[1]
    p_val <- param[2]
    n_val*p_val*(1 - p_val)*exp(tvec) / (1 - p_val + p_val*exp(tvec))^2
  },
  K3_vectorized_func = function(tvec, param) {
    n_val <- param[1]
    p_val <- param[2]
    tmp_ <- 1 - p_val + p_val * exp(tvec)
    n_val*(1 - p_val)*p_val*exp(tvec)*(1 - p_val - p_val * exp(tvec)) / tmp_^3
  },
  K4_vectorized_func = function(tvec, param) {
    n_val <- param[1]
    p_val <- param[2]
    tmp0 <- p_val*exp(tvec)
    tmp1 <- 1 - p_val
    tmp2 <- tmp1 + tmp0
    n_val*tmp1*tmp0*(4*p_val*tmp0 + tmp0^2 - 4*tmp0 + tmp1^2) / (tmp2^4)
  },
  
  
  analytic_tvec_hat_func = function(x, param) {
    n_val <- param[1]
    p_val <- param[2]
    
    # numeric check omitted for brevity, but typically you'd ensure y<n
    log(x*(1 - p_val)) - log(p_val*(n_val - x))
  },
  
  ############################
  op_name = "BinomialCGF"
)















#' @noRd
validateBinomialLengths <- function(vec, nVals, pVals, iidReps) {
  len_vec <- length(vec)
  ln <- length(nVals)
  lp <- length(pVals)
  
  if (ln != lp) stop(sprintf("Lengths of Binomial parameters nVals (%d) and pVals (%d) must match.", ln, lp))
  
  if (!is.null(iidReps)) {
    expected_len <- ln * iidReps
    if (len_vec != expected_len) {
      stop(sprintf(
        "Length of (tvec/x) is %d; expected %d (n/p length=%d, iidReps=%d).",
        len_vec, expected_len, ln, iidReps
      ))
    }
  } else if (len_vec %% ln != 0) {
    stop(sprintf(
      "Length mismatch: length(tvec or x)=%d, length(nVals)=%d. Either set 'iidReps' or ensure lengths match or are multiples.",
      len_vec, ln
    ))
  }
}


#' @noRd
.BinomialModelCGF_internal <- function(iidReps, ...){
  createCGF_fromVectorisedFunctions(
    
    K_vectorized_func = function(tvec, np) {
      # np => c( n_1,..., n_ln, p_1,..., p_ln ) from param_adaptor
      # first half => n, second half => p
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]

      validateBinomialLengths(tvec, nVals, pVals, iidReps)
      nVals*log(1 - pVals + pVals*exp(tvec)) 
    },
    
    
    
    K1_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      validateBinomialLengths(tvec, nVals, pVals, iidReps)
      nVals*pVals*exp(tvec) / (1 - pVals + pVals*exp(tvec))
    },
    
    
    
    K2_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      validateBinomialLengths(tvec, nVals, pVals, iidReps)
      
      nVals*pVals*(1 - pVals)*exp(tvec) / (1 - pVals + pVals*exp(tvec))^2
    },
    
    
    
    K3_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      validateBinomialLengths(tvec, nVals, pVals, iidReps)
      
      tmp_ <- 1 - pVals + pVals * exp(tvec)
      nVals*(1 - pVals)*pVals*exp(tvec)*(1 - pVals - pVals * exp(tvec)) / tmp_^3
    },
    
    
    
    K4_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      validateBinomialLengths(tvec, nVals, pVals, iidReps)
      
      tmp0 <- pVals*exp(tvec)
      tmp1 <- 1 - pVals
      tmp2 <- tmp1 + tmp0
      nVals*tmp1*tmp0*(4*pVals*tmp0 + tmp0^2 - 4*tmp0 + tmp1^2) / (tmp2^4)
    },
    
    
    
    ####################
    ## analytic_tvec_hat_func
    ####################
    analytic_tvec_hat_func = function(x, np) {
      # Solve K'(t_i)= y_i => t_i
      # Typically: t = log( y(1-p)/ [p(n-y)] )
      # but for a vector scenario, we do elementwise
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      validateBinomialLengths(x, nVals, pVals, iidReps)
      
      # elementwise
      # t_i = log( x_i(1-p_i) ) - log( p_i( n_i - x_i ) )
      # must ensure x_i< n_i for meaningful. 
      log(x * (1 - pVals)) - log(pVals*(nVals - x))
    },
    
    op_name = "BinomialModelCGF",
    ...
  )
}












#' Create a Parametric Binomial CGF Object
#'
#' @description
#' This creates a CGF object for Binomial(\eqn{n, prob}), where \eqn{n(\theta)} and
#' \eqn{prob(\theta)} can each be a function (or adaptor) of the user parameter vector \eqn{\theta}.
#' It supports i.i.d. replication (via \code{iidReps}) and non-identical usage by letting
#' \eqn{n/p} return vectors.
#'
#' @param n A function or adaptor returning \eqn{n} (an integer) or a vector of integers.
#' @param prob A function or adaptor returning \eqn{prob} in (0,1) or a vector of probabilities.
#' @param iidReps Either \code{"any"} or a positive integer specifying how many
#'   i.i.d. blocks are expected. Defaults to \code{"any"}, meaning no restriction on the length of \code{tvec}.
#' @param ... Additional arguments passed to CGF creating functions.
#'
#'
#' @return A \eqn{CGF} object.
#'
#' @examples
#' # ex 1) i.i.d. scenario: n=10, p=0.3 => param= c(10,0.3)
#' # using direct functions:
#' n_fn <- function(th) th[1]
#' p_fn <- function(th) th[2]
#' my_cgf <- BinomialModelCGF(n_fn, p_fn, iidReps=1) # if iidReps="any", length(tvec) determines the number of i.i.d. replicates
#' my_cgf$K1(0, c(10,0.3))
#'
#' # ex 2) non-identical scenario: n= c(5, 10), p= c(0.2, 0.7)
#' #   => param => c(5,10, 0.2,0.7)
#' n_adapt <- function(th) th[1:2] # OR n_adapt <- adaptor(indices = 1:2)
#' p_adapt <- function(th) th[3:4] # OR p_adapt <- adaptor(indices = 3:4)
#' my_cgf2 <- BinomialModelCGF(n_adapt, p_adapt, iidReps="any") # default iidReps="any", if iidReps=m, then m-i.i.d. blocks are expected => length(tvec) must be a multiple of m
#' # e.g. tvec= c(0,0) => length=2 => we have 2 binomial distributions: (5,0.2) and (10,0.7)
#' my_cgf2$K1(c(0,0), c(5,10,0.2,0.7))
#' 
#'
#' @export
BinomialModelCGF <- function(n, prob, iidReps = "any", ...) {
  if (is.character(iidReps) && length(iidReps) == 1 && tolower(iidReps) == "any") iidReps <- NULL
  if (!is.null(iidReps)) {
    if (length(iidReps) != 1 || is.infinite(iidReps) || !is.numeric(iidReps) ||
        iidReps < 1 || iidReps != as.integer(iidReps) )  {
      stop("'iidReps' must be 'any' or a positive integer.")
    }
  }
  
  n_fn <- validate_function_or_adaptor(n)
  p_fn <- validate_function_or_adaptor(prob)
  
  base_cgf <- .BinomialModelCGF_internal(iidReps, ...)
  
  
  adaptCGF(
    cgf = base_cgf,
    adaptor = function(theta) c(n_fn(theta), p_fn(theta))
  )
  
}











