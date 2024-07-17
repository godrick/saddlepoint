
CustomCGF <- R6::R6Class(
  "CustomCGF",
  public = list(
    Kvectorized = NULL,
    K1vectorized = NULL,
    K2vectorized = NULL,
    K3vectorized = NULL,
    K4vectorized = NULL,
    ineq_vectorized = NULL,
    
    initialize = function(Kvectorized, K1vectorized, K2vectorized, K3vectorized = NULL, 
                          K4vectorized = NULL, ineq_vectorized = NULL) {
      validate_function <- function(f, name) { # Helper function to validate user-provided functions
        if (!is.function(f)) { stop(sprintf("%s must be a function.", name)) }
        args <- names(formals(f))
        if (length(args) != 2 || !identical(args, c("tvec", "parameter_vector"))) {stop(sprintf("%s must have exactly two arguments named 'tvec' and 'parameter_vector'.", name))}
      }
      # Validate and assign mandatory functions
      validate_function(Kvectorized, "Kvectorized")
      validate_function(K1vectorized, "K1vectorized")
      validate_function(K2vectorized, "K2vectorized")
      self$Kvectorized = Kvectorized
      self$K1vectorized = K1vectorized
      self$K2vectorized = K2vectorized
      

      # Assign or default and validate optional functions
      self$K3vectorized <- if (is.null(K3vectorized)) {
        function(tvec, parameter_vector) rep(0, length(tvec))
      } else {
        validate_function(K3vectorized, "K3vectorized")
        K3vectorized
      }
      self$K4vectorized <- if (is.null(K4vectorized)) {
        function(tvec, parameter_vector) rep(0, length(tvec))
      } else {
        validate_function(K4vectorized, "K4vectorized")
        K4vectorized
      }
      self$ineq_vectorized <- if (is.null(ineq_vectorized)) {
        function(tvec, parameter_vector) numeric(0)
      } else {
        validate_function(ineq_vectorized, "ineq_vectorized")
        ineq_vectorized
      }
    },
    
    
    make_customCGF = function() {
      createCGF(make_CustomVectorizedScalarCGF(self$Kvectorized, self$K1vectorized, self$K2vectorized, 
                                               self$K3vectorized, self$K4vectorized, self$ineq_vectorized))
    }
  )
)




#' Create an additional CGF object.
#' 
#' This function generates a CGF object for univariate distributions, with vectorization to handle iid data.
#' All functions must take two arguments: tvec (vector for independent variable values) 
#' and parameter_vector (vector of underlying parameters).
#' 
#' @param Kvectorized A function representing the vectorized version of K.
#' @param K1vectorized A function representing the vectorized version of K1 (first derivative of K).
#' @param K2vectorized A function representing the vectorized version of K2 (second derivative of K).
#' @param K3vectorized A function representing the vectorized third derivative of K, essential for discrepancy calculations. If omitted, it defaults to returning a zero vector the length of tvec. Without an implementation of `K3vectorized` function, discrepancy calculations should not be performed.
#' @param K4vectorized A function representing the vectorized fourth derivative of K. If omitted, it defaults to returning a zero vector the length of tvec. Similar to `K3vectorized`, `K4vectorized` is essential for discrepancy calculations, which should not be performed without an implementation of this function.
#' @param ineq_vectorized A function representing the inequality constraints for the domain of the CGF. The constraints defined in this function should be non-positive for valid values of tvec. If a CGF does not require any constraints, it defaults to returning an empty vector.
#' 
#' @return An object of class 'CGF'.
#' 
#' @examples
#' \dontrun{
#' # Define the vectorized functions of the CGF object for a gamma distribution
#' vectorized_k <- function(tvec, parameter_vector){
#'   shape <- parameter_vector[1]
#'   rate <- parameter_vector[2]
#'   -shape * log1p(-tvec / rate)
#' }
#' vectorized_k1 <- function(tvec, parameter_vector){
#'   shape <- parameter_vector[1]
#'   rate <- parameter_vector[2]
#'   shape / (rate - tvec)
#' }
#' vectorized_k2 <- function(tvec, parameter_vector){
#'   shape <- parameter_vector[1]
#'   rate <- parameter_vector[2]
#'   shape / (rate - tvec)^2
#' }
#' vectorized_k3 <- function(tvec, parameter_vector){
#'   shape <- parameter_vector[1]
#'   rate <- parameter_vector[2]
#'   2 * shape / (rate - tvec)^3
#' }
#' vectorized_k4 <- function(tvec, parameter_vector){
#'   shape <- parameter_vector[1]
#'   rate <- parameter_vector[2]
#'   6 * shape / (rate - tvec)^4
#' }
#' # ineq_constraint is such that: tvec < rate
#' # The value of the next function is constrained to be non-positive: tvec - rate < 0
#' vectorized_ineq_constrain_fn <- function(tvec, parameter_vector){
#'   rate <- parameter_vector[2]
#'   tvec - rate
#' }
#'
#' # Create the CGF object for the gamma distribution
#' my_gamma_cgf <- createCustomCGF(Kvectorized = vectorized_k, K1vectorized = vectorized_k1, 
#'                                 K2vectorized = vectorized_k2, K3vectorized = vectorized_k3,
#'                                 K4vectorized = vectorized_k4, ineq_vectorized = vectorized_ineq_constrain_fn)
#'
#' # Example usage: Evaluate the first derivative at tvec = 0, parameters shape = 2, rate = 0.3
#' # Compare this result to the main CGF object GammaCGF implementation 
#' my_gamma_cgf$K1(0, c(2, 0.3)) == GammaCGF$K1(0, c(2, 0.3))
#'}
#' @export
createCustomCGF <- function(Kvectorized, K1vectorized, K2vectorized, 
                            K3vectorized = NULL, K4vectorized = NULL, ineq_vectorized = NULL) {
  custom_cgf <- CustomCGF$new(Kvectorized, K1vectorized, K2vectorized, K3vectorized, K4vectorized, ineq_vectorized)
  custom_cgf$make_customCGF()
}




