# '.saddlepointExternal_env' is an environment used to store and manage CGF objects loaded from external C++ code.
# It's intended for internal use only.
.saddlepointExternal_env <- new.env()


#' @title Compile and load externally created 'CGF' objects.
#'
#' @description
#' This function compiles and loads 'CGF' objects and associated functions that are not predefined in the 'saddlepoint' package.
#' The C++ code should be written following the provided template.
#'
#' @param file A character string giving the directory and name of the file to be compiled and loaded.
#' @param cacheDir A character string giving the directory to be used for caching compiled code.
#' If not provided, the current working directory will be used.
#' @param ... Other arguments passed on to `Rcpp::sourceCpp()`.
#'
#' @return A list containing all functions exported using `// [[Rcpp::export]]`. Each object is associated with its name in the list.`CGF_with_AD` objects will be of class 'CGF'.
#'
#' @export
sourceCppInSaddlepointEnv <- function(file, cacheDir = getwd(), ...) {
  if (!dir.exists(cacheDir)) dir.create(cacheDir, recursive = TRUE)

  rm(list = ls(envir = .saddlepointExternal_env), envir = .saddlepointExternal_env)

  # Compile and load the C++ code into .saddlepointExternal_env environment
  Rcpp::sourceCpp(file = file, cacheDir = cacheDir, env = .saddlepointExternal_env, ...)

  ### Some other related ideas -
  ### This function should load all exported C++ functions into the specific environment, wrap the ones that end with 'CGF', and then return these 'CGF' objects.
  ### Identify the CGF functions
  ### cgf_functions = grep("CGF$", ls(.saddlepointExternal_env), value = TRUE)

  cgf_functions = ls(.saddlepointExternal_env)

  for (fun_name in cgf_functions) { # Loop over the identified functions and create CGF objects
    fun = get(fun_name, envir = .saddlepointExternal_env)
    args = formals(fun) # Get function arguments
    ### cat(paste0("Loaded CGF function: '", fun_name, "' arguments: " ))
    cat(paste0("Loaded function: '", fun_name, "' arguments: " ))
    print(names(args))
    ### assign(fun_name, function(...) createUserCGF(do.call(fun, list(...))), envir = .saddlepointExternal_env)
    assign(fun_name, function(...) CGF(do.call(fun, list(...))), envir = .saddlepointExternal_env)
  }

  mget(ls(.saddlepointExternal_env), envir = .saddlepointExternal_env)
}


