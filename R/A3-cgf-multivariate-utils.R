
# # Validate block_size if provided (may be NULL for some CGFs)
# # #' @noRd
# .check_block_size <- function(block_size) {
#   if (is.null(block_size)) return(invisible())
#   if (!is.numeric(block_size) || length(block_size) != 1 || !is.finite(block_size) ||
#       block_size < 1 || block_size != as.integer(block_size)) {
#     stop("'block_size' must be a positive integer or NULL.")
#   }
#   invisible()
# }

# Validate block_size if provided (may be NULL, integer, or function(param) -> integer)
#' @noRd
.check_block_size <- function(block_size) {
  if (is.null(block_size)) return(invisible())
  if (is.function(block_size)) {
    # Optional: be strict about signature
    if (length(formals(block_size)) != 1L)
      stop("'block_size' function must have exactly one argument, e.g. function(param) ...")
    return(invisible())
  }
  if (!is.numeric(block_size) || length(block_size) != 1L || !is.finite(block_size) ||
      block_size < 1 || block_size != as.integer(block_size)) {
    stop("'block_size' must be a positive integer, a function(param)->integer, or NULL.")
  }
  invisible()
}

# Evaluate block_size for this call (coerces function to integer)
###### I might have a problem here // param might be advector and maybe what i need is some form of
###### RTMB:::getValues if param is advector
#' @noRd
.block_size_value <- function(block_size, param) {
  if (is.null(block_size)) return(NULL)
  if (is.function(block_size)) {
    return(block_size(param))
  }
  block_size
}


# Resolve replication layout for a given N = length(tvec), optional block_size, and iidReps
# Returns c(d = block_size, B = number_of_blocks)
#' @noRd
.resolve_rep_layout <- function(N, block_size, iidReps) {

  if (is.null(block_size)) {
    if (identical(iidReps, "any")) {
      # No replication enforced: treat as a single block of size N
      return(c(d = as.integer(N), B = 1L))
    } else {
      m <- as.integer(iidReps)
      if (N %% m != 0L) {
        stop(sprintf("length(tvec)=%d is not divisible by iidReps=%d.", N, m))
      }
      return(c(d = as.integer(N %/% m), B = m))
    }
  } else {
    d <- as.integer(block_size)
    if (identical(iidReps, "any")) {
      if (N %% d != 0L) {
        stop(sprintf("length(tvec)=%d is not a multiple of block_size=%d.", N, d))
      }
      return(c(d = d, B = as.integer(N %/% d)))
    } else {
      m <- as.integer(iidReps)
      if (N != d * m) {
        stop(sprintf("length(tvec)=%d != block_size * iidReps = %d * %d.", N, d, m))
      }
      return(c(d = d, B = m))
    }
  }
}



