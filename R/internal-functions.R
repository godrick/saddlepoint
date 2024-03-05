# #' Check for "tag" attribute in an object and its sub-attributes
# #'
# #' This function checks if the attribute "tag" is present in an object or any of its sub-attributes.
# #' The check is performed recursively, so it also considers sub-attributes of sub-attributes, and so on.
# #'
# #' @param x An object to check for the "tag" attribute.
# #' @return TRUE if the "tag" attribute is found, otherwise FALSE.
# hasTag <- function(x) {
#   if (!is.null(attributes(x)$tag)) return(TRUE)
#
#   if (is.list(x)) { # check each element of x
#     for (i in x) {
#       # Now recursively check for each element of the list
#       if (hasTag(i)) return(TRUE)
#     }
#   }
#
#   if (!is.null(attributes(x))) { # check attributes of attributes
#     for (i in attributes(x)) {
#       if (hasTag(i)) return(TRUE)
#     }
#   }
#   FALSE
# }


#' @title Create CppAD::ADfun Object Pointer for K1
#' @description This internal function creates an CppAD::ADFun object pointer for K1
#' from where we can extract the K1 value and gradient with respect to both 'tvec' and 'theta'.
#'
#' @param tvec (A numeric vector.) In estimation, the length of 'tvec' should match the number of observations.
#' @param theta (A numeric vector.) The length of 'theta' should match the number of parameters of interest in the model.
#' @param ModelCGF A CGF object formulated as a function of 'theta' using the any of the adaptor creating functions or model-based CGFs.
#'
#' @return An 'externalptr' which is internally CppAD::ADFun object pointer for K1.
#'
#' @keywords internal
#' @noRd
#'
#' @seealso \code{\link{adaptor}}
#'
#' @examples
#' \dontrun{
#' # ...
#' }
make.ADFunK1 <- function(tvec, theta, ModelCGF){
  if (!is.numeric(tvec) || !is.numeric(theta)) stop("tvec and theta must be numeric vectors")
  if (!is(ModelCGF, "CGF")) stop("'ModelCGF' is not defined for ", class(ModelCGF))
  makeADFunK1(tvec = tvec, theta = theta, modelCGF = ModelCGF)
}


# @title Create CppAD::ADfun Object Pointer for saddlepoint negative log-likelihood
# @description This internal function creates a CppAD::ADFun object pointer for the saddlepoint negative log-likelihood.
# from where we can extract the negative log-likelihood value and gradient with respect to both 'tvec' and 'theta'.
# 
# @param tvec (A numeric vector.) In estimation, the length of 'tvec' should match the number of observations.
# @param theta (A numeric vector.) The length of 'theta' should match the number of parameters of interest in the model.
# @param ModelCGF A CGF object formulated as a function of 'theta' using any of the adaptor creating functions or model-based CGFs.
# 
#@return An 'externalptr' which is internally a CppAD::ADFun object pointer for the saddlepoint negative log-likelihood.
# 
# @keywords internal
# @noRd
# 
# @seealso \code{\link{adaptor}}
# 
# @examples
# \dontrun{
# # ...
# }
# make.ADFunNegll <- function(tvec, theta, ModelCGF){
#   if (!is.numeric(tvec) || !is.numeric(theta)) stop("tvec and theta must be numeric vectors")
#   if (!is(ModelCGF, "CGF")) stop("'ModelCGF' is not defined for ", class(ModelCGF))
#
#   if ( hasTag(ModelCGF) ) {
#     # Call external function
#     if (!exists("makeADFunNegllExternal", envir = .saddlepointExternal_env)) stop("makeADFunNegllExternal function not found. Use function 'sourceCppInSaddlepointEnv' to create your CGF")
#     fun <- get("makeADFunNegllExternal", envir = .saddlepointExternal_env)
#     return(fun(tvec = tvec, theta = theta, modelCGF = ModelCGF))
#   }
#   # Call internal function
#   makeADFunNegll(tvec = tvec, theta = theta, modelCGF = ModelCGF)
# }
# make.ADFunNegll <- function(tvec, theta, ModelCGF){
#   if (!is.numeric(tvec) || !is.numeric(theta)) stop("tvec and theta must be numeric vectors")
#   if (!is(ModelCGF, "CGF")) stop("'ModelCGF' is not defined for ", class(ModelCGF))
#   makeADFunNegll(tvec = tvec, theta = theta, modelCGF = ModelCGF)
# }

#' @title Create CppAD::ADfun Object Pointer for Inequality Constraints
#' @description This internal function creates a CppAD::ADFun object pointer for the saddlepoint-based inequalities/constraints on 'tvec' that determine the validity of the 'ModelCGF'.
#' From where we can extract the constraints and gradients with respect to both 'tvec' and 'theta'.
#'
#' @param tvec (A numeric vector.) In estimation, the length of 'tvec' should match the number of observations.
#' @param theta (A numeric vector.) The length of 'theta' should match the number of parameters of interest in the model.
#' @param ModelCGF A CGF object formulated as a function of 'theta' using any of the adaptor creating functions or model-based CGFs.
#'
#' @return An 'externalptr' which is internally a CppAD::ADFun object pointer for the saddlepoint-based inequalities/constraints.
#'
#' @keywords internal
#' @noRd
#'
#' @seealso \code{\link{adaptor}}
#'
#' @examples
#' \dontrun{
#' # ...
#' }
make.ADFunIneqConstraint <- function(tvec, theta, ModelCGF){
  if (!is.numeric(tvec) || !is.numeric(theta)) stop("tvec and theta must be numeric vectors")
  if (!is(ModelCGF, "CGF")) stop("'ModelCGF' is not defined for ", class(ModelCGF))
  makeADFunIneqConstraint(tvec = tvec, theta = theta, modelCGF = ModelCGF)
}


