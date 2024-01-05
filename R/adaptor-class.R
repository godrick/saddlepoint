#' @title adaptor
#'
#' @description
#' The adaptor class is an externalptr object that points to related C++ objects.
#'
#' @slot .Data An externalptr object.
#' @return An adaptor object.
#'
#'
#'@noRd
setClass("adaptor", slots = list(.Data = "externalptr"))

#' adaptor Constructor
#'
#' This function constructs an object of class "adaptor" using a valid externalptr object.
#'
#' @param ptr An externalptr object. It should not be an empty externalptr object.
#'
#' @return An object of class "adaptor".
#' @importFrom methods is new
#' @seealso \code{\link{CGF}}
#' @examples
#' \dontrun{
#'   # assuming ptr is a valid externalptr object
#'   adaptor <- adaptor(ptr)
#' }
#' @noRd
adaptor <- function(ptr) {
  if (!is(ptr, "externalptr")) stop("ptr must be of type externalptr")
  if (isTRUE( identical(ptr, new("externalptr"))) ) stop("ptr is an empty externalptr object")
  structure(.Data = ptr, class = "adaptor")
}
