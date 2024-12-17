#' @title Sum of 'n' i.i.d random variables.
#'
#' @description
#' Computes the CGF object for the sum of 'n' i.i.d random variables.
#'
#'
#' @usage sumOfiidCGF(cgf, n)
#'
#' @param cgf A CGF object for the random variable to be summed.
#' @param n An integer representing the number of i.i.d random variables to be summed.
#' @return A CGF object.
#' @export
sumOfiidCGF <- function(cgf, n){
  if (!is(cgf, "CGF")) stop("'cgf' is not defined for ", class(cgf))
  if(!is.numeric(n) || n != as.integer(n) || n <= 0 || length(n) != 1) stop("n must be a single positive integer")
  createCGF(make_SumOfIIDCGF(cgf$get_ptr(), n))
}

#' @title Sum of Independent CGFs
#' 
#' @description
#' Computes the CGF object for the sum of independent (identical or non-identical) random variables.
#'
#' @param ... CGF objects: Any number of CGF objects. 
#' Each argument passed through `...` must be an object of class 'CGF', each corresponding to an independent random variable.
#' 
#' @return A CGF object for the sum of independent random variables.
#' @export
sumOfIndependentCGF <- function(...) {
  cgfs <- list(...)
  if (!length(cgfs)) stop("At least one 'CGF' object must be provided.")
  if ( any(vapply(cgfs, function(x) !is(x, "CGF"), FALSE)) ) stop("All arguments must be of class 'CGF'." )
  createCGF( make_SumOfIndependentCGF( lapply(cgfs, function(x) x$get_ptr()) ) )
}

#' @title CGF of a linearly mapped random variable.
#'
#' @description
#' Computes the CGF object of a linearly transformed random variable.
#'
#'
#' @param cgf The CGF object for the random variable before the linear transformation.
#' @param matrix_A The transformation matrix applied to the random variable.
#' @return A CGF object
#' @export
linearlyMappedCGF <- function(cgf, matrix_A) {
  if (!is(cgf, "CGF"))  stop("'cgf' is not defined for ", class(cgf))
  if (!is.matrix(matrix_A)) stop("'matrix_A' is not defined for ", class(cgf))
  createCGF(make_LinearlyMappedCGF(cgf$get_ptr(), matrix_A))
}

#' @title CGF of a randomly stopped sum.
#' @description 
#' Computes the CGF object for a randomly stopped sum \eqn{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}{Y = sum(X[1], X[2], ..., X[N])}, where \eqn{N} and \eqn{\tilde{X}_i} are random variables. Here, \eqn{\tilde{X}_i} are i.i.d. random variables independent of \eqn{N}, which is a non-negative integer-valued random variable.
#' 
#' @details
#' The setup is as follows:
#' \deqn{{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}}
#' The random variable \eqn{\tilde{Y}} is computed as the sum of \eqn{N} i.i.d copies of the random vector \eqn{\tilde{X}_i}. Each \eqn{\tilde{X}_i} is a random vector of arbitrary dimension \eqn{d}.
#' For example, when \eqn{d=3}:
#' \deqn{ \tilde{X_1} = (X_{11}, X_{12}, X_{13}) }
#' \deqn{ \tilde{X_2} = (X_{21}, X_{22}, X_{23}) }
#' \deqn{ ... }
#' \deqn{ \tilde{X_N} = (X_{N1}, X_{N2}, X_{N3}) }
#' Here, \eqn{X_{i1}, X_{i2}, X_{i3}} are scalars but the vectors \eqn{\tilde{X_i}} are i.i.d copies.
#'
#' @param count_cgf The CGF object representing the random variable \eqn{N}.
#' @param summand_cgf The CGF object representing the random variables \eqn{\tilde{X}_i}.
#' @return A CGF object of the random variable \eqn{\tilde{Y}}.
#' @export
randomlyStoppedSumCGF <- function(count_cgf, summand_cgf) {
  if (!is(count_cgf, "CGF")) stop("'count_cgf' is not defined for ", class(count_cgf))
  if (!is(summand_cgf, "CGF")) stop("'summand_cgf' is not defined for ", class(summand_cgf))
  createCGF(make_RandomlyStoppedSumCGF(count_cgf = count_cgf$get_ptr(), summand_cgf = summand_cgf$get_ptr()))
}

#' Transforms a CGF object to an i.i.d-based CGF object.
#'
#'
#' @param cgf (CGF object): The CGF for the random vector of dimension (block_size).
#' @param block_size (positive integer): The vector length of a single replicate.
#' @return A CGF object.
#' @export
iidReplicatesCGF <- function(cgf, block_size) {
  if (!is(cgf, "CGF")) stop("'cgf' is not defined for ", class(cgf))
  if (!is.numeric(block_size) || block_size != as.integer(block_size) || block_size <= 0 || length(block_size) != 1) stop("'block_size' must be a positive integer")
  createCGF(make_IIDReplicatesCGF(cgf$get_ptr(), block_size))
}


#' @title Concatenation of CGF objects.
#'
#' @description Computes the CGF object for a concatenated vector of independent random variables, with each represented by its respective CGF object and a corresponding dimension. It facilitates the modeling of a concatenated random vector, where each segment of the vector is specified by its CGF object and a corresponding dimension.
#'
#' @param ... A sequence of alternating CGF objects and integers. Each CGF object must be immediately followed by its corresponding vector length, which denotes the dimension of the underlying random variable. 
#'
#' @return A CGF object.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
concatenationCGF <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]])) args <- args[[1]]
  if (length(args) == 0) stop("At least one CGF object and its corresponding length must be provided.")
  if (length(args) %% 2 != 0) stop("Arguments must be provided in pairs of CGF objects and their lengths.")
  if (length(args) / 2 < 2) stop("At least two pairs of CGF objects and their lengths must be provided to perform concatenation.")
  cgf_ptrs_and_lengths <- vector("list", length(args) / 2)
  for (i in seq(1, length(args), by = 2)) {
    cgf <- args[[i]]
    length_ <- args[[i + 1]]
    if (!inherits(cgf, "CGF")) stop("Argument at position ", i, " must be a CGF object.")
    if (!is.numeric(length_) || length_ != as.integer(length_) || length_ <= 0) stop("Argument at position ", i + 1, " must be a single positive integer.")
    
    cgf_ptrs_and_lengths[[i / 2 + 0.5]] <- list(cgf$get_ptr(), length_)
  }
  createCGF(make_ConcatenationCGF(cgf_ptrs_and_lengths))
  # if (!is.list(cgf_with_lengths)) stop("'cgf_with_lengths' must be a list, got ", class(cgf_with_lengths), ".")
  # if (length(cgf_with_lengths) <= 1) stop("'cgf_with_lengths' must contain more than one pair of CGF and length.")
  # cgf_ptrs_and_lengths <- lapply(cgf_with_lengths, function(item) {
  #   if (!is.list(item) || length(item) != 2)  stop("Each pair/sublist must be a list of exactly two elements (CGF and length).")
  #   if (!is(item[[1]], "CGF")) stop("The first element of each pair/sublist must be a CGF object.")
  #   if (!is.numeric(item[[2]]) || item[[2]] != as.integer(item[[2]]) || item[[2]] <= 0 || length(item[[2]]) != 1) stop("The second element of each pair/sublist in 'cgfWithLengths' must be a single positive integer.")
  #   list(ptr = item[[1]]$get_ptr(), length = item[[2]])
  # })
  # createCGF(make_ConcatenationCGF(cgf_ptrs_and_lengths))
}

### Quick check for the ConcatenationCGF function
### TO DO: Rewrite the example to use the new functions/objects
# if(FALSE){
#   cgf1 = PoissonModelCGF(lambda = adaptorUsingIndices(1))
#   cgf2 = BinomialModelCGF(n = adaptorUsingIndices(2), prob = adaptorUsingIndices(3))
#   concat_cgf = ConcatenationCGF(cgfWithLengths = list(list(cgf1,3), list(cgf2,2)))
#
#
#   A0 = rbind(diag(3), matrix(0, 2, 3))
#   A1 = rbind(matrix(0, 3, 2), diag(2))
#   cgf1_alt = linearlyMappedCGF(PoissonCGF, A0)
#   cgf2_alt = linearlyMappedCGF(BinomialCGF, A1)
#   cgf1_alt1 = adaptCGF(cgf1_alt, adaptorUsingIndices(1))
#   cgf2_alt2 = adaptCGF(cgf2_alt, adaptorUsingIndices(2:3))
#   my.cgf = sumOfIndependentCGF(cgf1_alt1, cgf2_alt2)
#
#
#   K(c(0.1, 0.2, 0.3, 0.4, 0.5), c(12, 10, 0.64), concat_cgf) == K(c(0.1, 0.2, 0.3, 0.4, 0.5), c(12, 10, 0.64), my.cgf)
# }











