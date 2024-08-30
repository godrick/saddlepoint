# A version of the RTMB internal code RTMB:::xtra
# We avoid using the RTMB:::xtra environment and instead use the exported RTMB::ADoverload function.
# This environment overrides the R primitives ([<-, diag<-, and c) with custom functions that internally use 
# RTMB::ADoverload to dynamically apply RTMB::advector-based operators.
# Usage: 
# 1. We can simply set rtmb_Xtra as the environment of a function that requires the custom definitions. (##### check this)
#     e.g. f1 <- function() {
#   		     # some operation using "[<-", "diag<-", and "c" 
#          }
#          environment(f) <- myXtra
#    ##### But this might lead to complications if some internal functions are in a different environment.
# 2. temporary usage: attach rtmb_Xtra during the execution of a function and detach it on exit. e.g.,
#    f2 <- function(x) {
#	    attach(myXtra, name = "myXtra", pos = 2)
#	    on.exit(detach("myXtra"), add = TRUE)
#	    # operations that require the custom definitions
#    }
# There maybe other options ##### check if "with" can be used.
## We use the environment in adaptedCGFs.R file when a user supplies an R function as an adaptor.
## This should also be used in cgf_from_rfunction.R file.
rtmb_Xtra <- local({
  "[<-" <- function(x, ..., value) {
    RTMB::ADoverload("[<-") (x, ..., value=value)
  }
  "diag<-" <- function(x, value) {
    RTMB::ADoverload("diag<-") (x, value)
  }
  c <- function(...) {
    RTMB::ADoverload("c")(...)
  }
  environment()
})
