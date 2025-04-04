#' @keywords package
#' @section Package development: 
#' The package is still in development and there may still be bugs and errors. 
#' While we do not expect major changes to the general user interface, there may be changes to the internal workings of the package as well as new additions and functionality.
"_PACKAGE"

## usethis namespace: start
#' @import RTMB
#' @import Matrix
#' @useDynLib saddlepoint, .registration = TRUE
## usethis namespace: end
NULL

# Additionally, the package effectively handles constraints in optimization problems and offers integration with C++ for accurate and efficient gradient computations. Its design allows for extensibility, accommodating user-defined 
# CGFs and operations.    This package is particularly useful for statisticians and data scientists who 
# work with random variables whose probability distribution functions are 
# computationally intractable. It streamlines the process of handling these 
# variables, making the analysis more accessible and manageable.