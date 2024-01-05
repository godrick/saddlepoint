#' @title CGF
#'
#' @description
#' Refer to the CGF() function for details.
#'
#' @slot .Data An externalptr object.
#'
#'
#'@noRd
setClass("CGF", slots = list(.Data = "externalptr"))

#' CGF Constructor
#'
#' This function constructs an object of class "CGF" using a valid externalptr object.
#'
#' @param ptr An externalptr object. It should not be an empty externalptr object.
#'
#' @return An object of class "CGF".
#' @examples
#' \dontrun{
#'   # assuming ptr is a valid externalptr object
#'   cgf <- CGF(ptr)
#' }
#' @importFrom methods is new
#' @export
CGF <- function(ptr) {
  if (!is(ptr, "externalptr")) stop("'ptr' must be of type externalptr")
  if (isTRUE( identical(ptr, new("externalptr"))) ) stop("ptr is an empty externalptr object")
  structure(.Data = ptr, class = "CGF")
}



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

