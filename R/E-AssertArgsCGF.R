#' @noRd
assertArgsCGF <- function(cgf, assertArgs_fn){
  
  # default checker will only try to validate iidReps
  # So it should only be checking that length of tvec is divisible by iidReps
  # so this will not check anything about the parameters, if that is needed, then we'd need to 
  # override the default checker
  checker_fn <- function(cgf, tvec, iidReps){
    if (length(tvec) %% iidReps != 0) stop("Length of tvec must be divisible by iidReps")
    
  }
  
  
}