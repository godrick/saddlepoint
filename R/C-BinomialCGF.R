# R/BinomialCGF.R






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
  
  
  analytic_tvec_hat_func = function(y, param) {
    n_val <- param[1]
    p_val <- param[2]
    
    # numeric check omitted for brevity, but typically you'd ensure y<n
    log(y*(1 - p_val)) - log(p_val*(n_val - y))
  },
  
  ############################
  op_name = "BinomialCGF"
)














#' @noRd
expandNandP <- function(vec, nVals, pVals, iidReps) {
  # Suppose 'vec' is tvec or y
  len_vec <- length(vec)
  ln <- length(nVals)
  lp <- length(pVals)
  
  if (ln != lp) stop(sprintf("Lengths of nVals (%d) and pVals (%d) must match.", ln, lp))
  
  
  # If we have an integer iidReps > 0
  if (!is.null(iidReps)) {
    expected_len <- ln * iidReps
    if (len_vec != expected_len) {
      stop(sprintf(
        "Length of 'tvec/observations' is %d; expected %d (n/p length=%d, iidReps=%d).",
        len_vec, expected_len, ln, iidReps
      ))
    }
    # # replicate each n, p across blocks
    # # We avoid lists for now
    # return(list(
    #   n_expanded = rep(nVals, times = iidReps),
    #   p_expanded = rep(pVals, times = iidReps)
    # ))
    return(
      c(rep(nVals, times = iidReps), rep(pVals, times = iidReps))
    )
  }
  
  # If iidReps is NULL, we allow more flexible usage
  if ((is.null(iidReps)) && (len_vec %% ln == 0) ) {
    # replicate if needed
    times_ <- len_vec / ln
    return(
      c(rep(nVals, times = times_), rep(pVals, times = times_))
    )
  }
  
  stop(sprintf(
    "Length mismatch: length(tvec/observations)=%d, lengths(nVals/pVals)=%d. Either set 'iidReps', or ensure lengths match or are factors.",
    len_vec, ln
  ))
}















#' Create a Parametric Binomial CGF Object
#'
#' @description
#' This creates a CGF object for Binomial(\eqn{n, prob}), where \eqn{n(\theta)} and
#' \eqn{prob(\theta)} can each be a function (or adaptor) of the user parameter vector \eqn{\theta}.
#' It supports i.i.d. replication (via \code{iidReps}) and non-identical usage by letting
#' n/p return vectors.
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
#' my_cgf <- BinomialModelCGF(n_fn, p_fn, iidReps=1) # if iidReps=NULL, length(tvec) determines the number of i.i.d. replicates
#' my_cgf$K1(0, c(10,0.3))
#'
#' # ex 2) non-identical scenario: n= c(5,10), p= c(0.2,0.7)
#' #   => param => c(5,10, 0.2,0.7)
#' n_adapt <- function(th) th[1:2] # OR n_adapt <- adaptor(indices = 1:2)
#' p_adapt <- function(th) th[3:4] # OR p_adapt <- adaptor(indices = 3:4)
#' my_cgf2 <- BinomialModelCGF(n_adapt, p_adapt, iidReps=NULL) # default iidReps=NULL, if iidReps=m, then m-i.i.d. blocks are expected => length(tvec) must be a multiple of m
#' # e.g. tvec= c(0,0) => length=2 => we have 2 binomial distributions: (5,0.2) and (10,0.7)
#' my_cgf2$K1(c(0,0), c(5,10,0.2,0.7))
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
  
  
  base_cgf <- createCGF_fromVectorisedFunctions(
    
    
    
    K_vectorized_func = function(tvec, np) {
      # np => c( n_1,..., n_ln, p_1,..., p_ln ) from param_adaptor
      # first half => n, second half => p
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(tvec, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      n_expanded*log(1 - p_expanded + p_expanded*exp(tvec)) 
    },
    

    
    K1_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(tvec, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      n_expanded*p_expanded*exp(tvec) / (1 - p_expanded + p_expanded*exp(tvec))
    },
    

    
    K2_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(tvec, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      n_expanded*p_expanded*(1 - p_expanded)*exp(tvec) / (1 - p_expanded + p_expanded*exp(tvec))^2
    },
    
    
    
    K3_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(tvec, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      tmp_ <- 1 - p_expanded + p_expanded * exp(tvec)
      n_expanded*(1 - p_expanded)*p_expanded*exp(tvec)*(1 - p_expanded - p_expanded * exp(tvec)) / tmp_^3
    },
    
    
    
    K4_vectorized_func = function(tvec, np) {
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(tvec, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      tmp0 <- p_expanded*exp(tvec)
      tmp1 <- 1 - p_expanded
      tmp2 <- tmp1 + tmp0
      n_expanded*tmp1*tmp0*(4*p_expanded*tmp0 + tmp0^2 - 4*tmp0 + tmp1^2) / (tmp2^4)
    },
    
    
    
    ####################
    ## analytic_tvec_hat_func
    ####################
    analytic_tvec_hat_func = function(y, np) {
      # Solve K'(t_i)= y_i => t_i
      # Typically: t = log( y(1-p)/ [p(n-y)] )
      # but for a vector scenario, we do elementwise
      ln <- length(np)/2
      idx_n <- 1:ln
      idx_p <- ln + idx_n
      nVals <- np[idx_n]
      pVals <- np[idx_p]
      
      out_ <- expandNandP(y, nVals, pVals, iidReps)
      n_expanded <- out_[idx_n]
      p_expanded <- out_[idx_p]
      
      # elementwise
      # t_i = log( y_i(1-p_i) ) - log( p_i( n_i - y_i ) )
      # must ensure y_i< n_i for meaningful. 
      log(y * (1 - p_expanded)) - log(p_expanded*(n_expanded - y))
    },
    
    op_name = "BinomialModelCGF",
    ...
  )
  
  # Finally, we adapt param => c(n(theta), p(theta)), i.e. the user-supplied n, p
  # => a single vector c( nvals, pvals ).
  # each is possibly a vector if the user wants non-identical.
  adaptCGF(
    cgf = base_cgf,
    param_adaptor = function(theta) {
      # Evaluate the user functions
      # or adaptors => numeric vectors
      c(n_fn(theta), p_fn(theta))
    }
  )
  
}











