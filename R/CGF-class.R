#' @importFrom methods is new
CGF <- R6::R6Class("CGF",
                   private = list(
                     ptr = NULL  # Holds the externalptr
                   ),
                   public = list(
                     initialize = function(ptr) {
                       if (!is(ptr, "externalptr"))  stop("ptr must be of type externalptr", call. = FALSE)
                       if (isTRUE(identical(ptr, new("externalptr")))) stop("ptr is an empty externalptr object")
                       private$ptr <- ptr
                     },
                     get_ptr = function() { private$ptr }, # Expose the externalptr directly
                     # print = function(...) {
                     #   cls <- class(self)
                     #   cls <- cls[cls != "R6"]
                     #   cat("An object of class '", paste(cls, collapse = " "), "'\n", sep = "")
                     # },
                     
                     
                     K = function(tvec, parameter_vector) {K_impl(tvec, parameter_vector, private$ptr)},
                     K1 = function(tvec, parameter_vector) {K1_impl(tvec, parameter_vector, private$ptr)},
                     K2 = function(tvec, parameter_vector) {K2_impl(tvec, parameter_vector, private$ptr)},
                     ineq_constraint = function(tvec, parameter_vector)  {ineq_constraint_impl(tvec, parameter_vector, private$ptr)}
                     # add: hat_f, neg_ll
                     
                     
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