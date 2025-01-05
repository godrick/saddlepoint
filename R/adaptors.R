# R/adaptors.R
# Objects: 

#### Check if we should retain this adaptor class
#### If we do then for each "adaptor" object, the call evaluateAdaptor(adaptor_obj, theta) should be enough! 


# #' @title adaptor
# #'
# #' @description
# #' An S4 class used to encapsulate how a parameter is extracted or computed
# #' from a parameter vector \code{theta}. This can be:
# #' \itemize{
# #'   \item A subset of \code{theta} by index,
# #'   \item A fixed numeric value,
# #'   \item A user-defined function.
# #' }
# #'
# #' Internally, the `.Data` slot is a plain R function of the form
# #' \code{function(param_vector) -> numeric}.
# #'
# #' @slot .Data A function(\code{param_vector}) -> numeric, implementing the logic.
# #' @noRd
# setClass("adaptor", slots = list(.Data = "function"))





# #' Create an adaptor object
# #'
# #' This function constructs an \code{adaptor} object from exactly one of:
# #' \enumerate{
# #'   \item \code{indices}: A numeric vector of positive integer indices, used to subset
# #'     \code{theta}.
# #'   \item \code{fixed_param}: A numeric value (or vector) to be returned regardless of \code{theta}.
# #'   \item \code{r_func}: A custom R function of the form
# #'     \code{function(theta) -> numeric}.
# #' }
# #'
# #' @param indices A numeric vector of positive integers.
# #' @param fixed_param A numeric vector or scalar to return.
# #' @param r_func A function(\code{theta}) -> numeric.
# #'
# #' @details
# #' \strong{Exactly one} of these three arguments must be non-\code{NULL}.
# #'
# #'
# #' @return An object of class \code{adaptor} whose \code{.Data} slot is a function.
# #'
# #'
# #' @export
# adaptor <- function(indices = NULL, fixed_param = NULL, r_func = NULL) {
#   args <- as.list(sys.call())[-1] # remove the function name from the captured call
#   if (length(args) != 1) stop("Please specify exactly one of 'indices', 'fixed_param', or 'r_func'")
#   if (is.null(names(args))) stop("Please specify the argument by its name")
#   
#   
#   if (names(args) == "indices") {
#     if (!is.numeric(indices) || any(indices != as.integer(indices)) || any(indices <= 0)) stop("'indices' must be positive integers")
#     func_ <- function(param_vector) {
#       param_vector[indices]
#     }
#     return(methods::new("adaptor", .Data = func_))
#   }
#   
#   if (names(args) == "fixed_param"){
#     if(!is.numeric(fixed_param)) stop(" 'fixed_param' is not defined for ", class(fixed_param)) 
#     func_ <- function(param_vector) fixed_param
#     return(methods::new("adaptor", .Data = func_))
#   }
#   
#   if (names(args) == "r_func"){
#     if(!is.function(r_func)) stop(" 'r_func' is not defined for ", class(r_func)) 
#     if (length(formals(r_func)) != 1) stop("The function must have exactly one argument.")
#     return(methods::new("adaptor", .Data = r_func))
#   }
#   
#   stop("Internal error: no valid argument chosen.")
# }



# #' Evaluate an \code{adaptor} object with a parameter vector
# #'
# #' @param adaptor An object of class \code{adaptor}.
# #' @param param_vector A numeric vector of parameters.
# #' @return The numeric result of calling \code{ad@.Data(param_vector)}.
# #' @noRd
# evaluateAdaptor <- function(adaptor, param_vector) {
#   if (!inherits(adaptor, "adaptor")) {
#     stop("Not a valid 'adaptor' object.")
#   }
#   fn <- adaptor@.Data
#   fn(param_vector)
# }






