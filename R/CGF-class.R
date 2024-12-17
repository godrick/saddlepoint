#' A CGF object
#' 
#' An object with methods for the cumulant generating function (CGF) and its related functions.
#' @details
#' A CGF object is a reference to an external C++ object that contains the cumulant generating function and its first- and second-order derivatives for a specific distribution. This object also supports functions for higher-order derivatives and additional features such as an inequality constraint function to manage domain constraints on the CGF. 
#' 
#' @references in preparation
#' 
#' @importFrom methods is new
#' @noRd
CGF <- R6::R6Class("CGF",
                   private = list(
                     ptr = NULL  # Holds the externalptr
                   ),
                   public = list(
                     #' @description 
                     #' Initializes a new CGF object with a valid external pointer.
                     #' @param ptr An external pointer 'externalptr' expected by the CGF class.
                     #' @return An object of class 'CGF'.
                     initialize = function(ptr) {
                       if (!is(ptr, "externalptr"))  stop("ptr must be of type externalptr", call. = FALSE)
                       if (isTRUE(identical(ptr, new("externalptr")))) stop("ptr is an empty externalptr object")
                       private$ptr <- ptr
                     },
                     #' @description
                     #' Returns the external pointer 'externalptr' associated with the CGF object.
                     #' @return An external pointer 'externalptr'.
                     get_ptr = function() { private$ptr }, 
                                   # print = function(...) {
                                   #   cls <- class(self)
                                   #   cls <- cls[cls != "R6"]
                                   #   cat("An object of class '", paste(cls, collapse = " "), "'\n", sep = "")
                                   # },
                     
                     #' @description The cumulant generating function 
                     #' @param tvec A numeric vector 
                     #' @param parameter_vector A numeric vector
                     #' @return A scalar value
                     K = function(tvec, parameter_vector) {K_impl(tvec, parameter_vector, private$ptr)},
                     #' @description The first derivative of the cumulant generating function
                     #' @param tvec A numeric vector
                     #' @param parameter_vector A numeric vector
                     #' @return A numeric vector
                     K1 = function(tvec, parameter_vector) {K1_impl(tvec, parameter_vector, private$ptr)},
                     #' @description The second derivative of the cumulant generating function
                     #' @param tvec A numeric vector
                     #' @param parameter_vector A numeric vector
                     #' @return A matrix
                     K2 = function(tvec, parameter_vector) {K2_impl(tvec, parameter_vector, private$ptr)},
                     #' @description Constraint function for the domain of the CGF
                     #' @param tvec A numeric vector
                     #' @param parameter_vector A numeric vector
                     #' @return A numeric vector
                     ineq_constraint = function(tvec, parameter_vector)  {ineq_constraint_impl(tvec, parameter_vector, private$ptr)},
                     # add: hat_f, neg_ll
                     saddlepoint_density = function(y, parameter_vector)  {
                       stop("saddlepoint_density not implemented yet")
                       # #### TO DO: check if this function can be implemented in C++
                       # saddlepoint_vals <- saddlepoint.eqn.solve(theta = parameter_vector, y = y, cgf = self,
                       #                                           sadd.eqn.opts = list(ftol_abs = 0, maxeval = 1e3, xtol_rel = 1.0e-12, print_level = 0) )
                       ## IF (multivariate) exp(-neg_ll_impl(saddlepoint_vals, parameter_vector, private$ptr))
                       # exp(-sapply(saddlepoint_vals, neg_ll_impl, parameter_vector, private$ptr))
                     }
                     
                     
                   )
)

#' Create a CGF object
#'
#' This function initializes a new CGF object with a valid external pointer.
#'
#' @param ptr An external pointer 'externalptr' expected by the CGF class.
#'
#' @return An object of class 'CGF'.
#' @noRd
createCGF <- function(ptr) {
  CGF$new(ptr)
}

# #' @title CGF
# #'
# #' @description
# #' Refer to the CGF() function for details.
# #'
# #' @slot .Data An externalptr object.
# #'
# #'
# #'@noRd
# setClass("CGF", slots = list(.Data = "externalptr"))



# CGF <- function(ptr) {
#   if (!is(ptr, "externalptr")) {
#     stop("'ptr' must be of type externalptr")
#   }
#   if (isTRUE(identical(ptr, new("externalptr")))) {
#     stop("ptr is an empty externalptr object")
#   }
#   obj <- structure(.Data = ptr, class = "CGF")
#   
#   
#   K <- function(tvec, parameter_vector) {
#     K_impl(tvec, parameter_vector, obj@.Data)
#   }
#   
#   attr(obj, "K") <- K
#   obj
# }


# CGF <- function(ptr) {
#   if (!is(ptr, "externalptr")) stop("'ptr' must be of type externalptr")
#   if (isTRUE( identical(ptr, new("externalptr"))) ) stop("ptr is an empty externalptr object")
#   structure(.Data = ptr, class = "CGF")
# }



# #' @title External CGF Constructor
# #'
# #' @description
# #' This function constructs an object of class "CGF" using a valid externalptr object. This function is used within `sourceCppInSaddlepointEnv` function to wrap externally created CGF objects.
# #'
# #'
# #' @param ptr An externalptr object. It should not be an empty externalptr object. TO DO: document more on this
# #'
# #' @return An object of class "CGF".
# #' @examples
# #' \dontrun{
# #'   # assuming ptr is a valid externalptr object
# #'   cgf <- createUserCGF(ptr)
# #' }
# #' @keywords internal
# createUserCGF <- function(ptr) {
#   cgf <- CGF(ptr)
#   attr(cgf, "tag") = "UserCGF"
#   cgf
# }

# Rcpp::loadModule("thetaGradientFunctions", TRUE)