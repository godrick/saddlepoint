#' @title CGF of a sum of 'n' i.i.d random variables.
#'
#' @description
#' Computes the CGF object for the sum of 'n' i.i.d random variables.
#'
#'
#' @usage sumOfiidCGF(baseCGF, n)
#'
#' @param baseCGF A CGF object the single random variable.
#' @param n An integer representing the number of i.i.d random variables to be summed.
#' @return A CGF object representing the sum of 'n' i.i.d random variables.
#' @export
sumOfiidCGF <- function(baseCGF, n){
  if (!is(baseCGF, "CGF")) stop("'baseCGF' is not defined for ", class(baseCGF))
  if(!is.numeric(n) || n != as.integer(n) || n <= 0 || length(n) != 1) stop("n must be a single positive integer")
  CGF(make_sum_of_iidCGF(baseCGF, n))
}
#' @title CGF of a sum of independent random variables.
#'
#' @description
#' Evaluates the CGF object for the sum of independent(identical or non-identical) random variables.
#'
#'
#' @usage sumOfIndependentCGF(baseCGF_list)
#'
#' @param baseCGF_list List of an arbitrary number of CGF objects, each corresponding to an independent random variable.
#' @return A CGF object representing the sum of independent random variables.
#' @export
sumOfIndependentCGF <- function(baseCGF_list) {
  if(!is.list(baseCGF_list)) stop(" 'baseCGF_list' is not defined for ", class(baseCGF_list))
  if (!length(baseCGF_list))  stop("At least one 'CGF' object must be provided.")
  if ( any(vapply(baseCGF_list, function(x) !is(x, "CGF"), FALSE)) ) stop("All input arguments must be of class 'CGF'." )
  CGF(make_sum_of_independentCGF(baseCGF_list))
}
#' @title CGF of a linearly mapped random variable.
#'
#' @description
#' Evaluates the CGF object of a linearly transformed random variable.
#'
#'
#' @param baseCGF The CGF object representing the random variable before the linear transformation.
#' @param matrix_A The transformation matrix applied to the random variable.
#' @return A CGF object representing the transformed random variable.
#' @export
linearlyMappedCGF <- function(baseCGF, matrix_A) {
  if (!is(baseCGF, "CGF")) stop("baseCGF must be of class 'CGF'.")
  if (!is.matrix(matrix_A)) stop("matrix_A must be a matrix.")
  CGF(make_linearly_mappedCGF(baseCGF, matrix_A))
}



#' @title CGF of a randomly stopped sum.
#' @description Evaluates the CGF for a randomly stopped sum \eqn{\tilde{Y} = \sum_{i=1}^{N} \tilde{X}_i}{Y = sum(X[1], X[2], ..., X[N])}, where \eqn{N} and \eqn{\tilde{X}_i} are random variables.
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
#' @param countCGF The CGF object representing the random variable \eqn{N}.
#' @param summandCGF The CGF object representing the random variables \eqn{\tilde{X}_i}.
#' @return A CGF object of the random variable \eqn{\tilde{Y}}.
#' @export
randomlyStoppedSumCGF <- function(countCGF, summandCGF) {
  if (!is(countCGF, "CGF")) stop("'countCGF' must be of class 'CGF'.")
  if (!is(summandCGF, "CGF")) stop("'summandCGF' must be of class 'CGF'.")
  CGF(make_randomly_stopped_sumCGF(count_base_cgf = countCGF, summand_base_cgf = summandCGF))
}

#' Transforms a scalar CGF to an i.i.d CGF
#'
#' @param scalarCGF A CGF object representing the scalar CGF to be transformed.
#' @return A CGF object representing the resulting i.i.d CGF.
#' @export
scalarToiidCGF <- function(scalarCGF) {
  if (!is(scalarCGF, "CGF")) stop("'scalarCGF' must be of class 'CGF'.")
  CGF(make_scalar_to_iidCGF(scalarCGF))
}

#' Transforms a CGF object to an i.i.d-based CGF object.
#'
#'
#' @param baseCGF (CGF object): The CGF for the random vector of dimension (block_size).
#' @param block_size (positive integer): The vector length of a single replicate.
#' @return A CGF object.
#' @export
iidReplicatesCGF <- function(baseCGF, block_size) {
  if (!is(baseCGF, "CGF")) stop("baseCGF must be of class 'CGF'")
  if (!is.numeric(block_size) || block_size != as.integer(block_size) || block_size <= 0 || length(block_size) != 1) stop("block_size must be a positive integer")
  CGF(make_iidReplicatesCGF(base_cgf = baseCGF, block_size = block_size))
}


#' @title Concatenation CGF
#' @description Enables the combination of multiple CGFs into a concatenated structure where each CGF has an associated vector length.
#'
#' @param cgfWithLengths A list where each entry consists of a CGF object and its associated vector length. Each element is a sublist containing:
#'   - The CGF object (CGF object).
#'   - The associated vector length for that CGF (positive integer).
#'
#' @return A CGF object.
#'
#' @examples
#' \dontrun{
#' # TO DO
#' }
#'
#' @export
concatenationCGF <- function(cgfWithLengths) {
  if (!is.list(cgfWithLengths)) stop("'cgfWithLengths' is not defined for class ", class(cgfWithLengths))
  if (length(cgfWithLengths) <= 1) stop("'cgfWithLengths' must be a list of more than one pair of CGF and length.")
  for (item in cgfWithLengths) {
    if (!is.list(item) || length(item) != 2)  stop("Each element in 'cgfWithLengths' must be a list of length 2.") # Each item is a list of length 2
    if (!is(item[[1]], "CGF")) stop("The first element of each sublist in 'cgfWithLengths' must be of class 'CGF'.") # Check the first element of each sublist
    if (!is.numeric(item[[2]]) || item[[2]] != as.integer(item[[2]]) || item[[2]] <= 0 || length(item[[2]]) != 1) stop("The second element of each sublist in 'cgfWithLengths' must be a positive integer.")
  }
  CGF(make_ConcatenationCGF(cgfWithLengths))
}
### Quick check for the ConcatenationCGF function
# if(F){
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
#   my.cgf = sumOfIndependentCGF(baseCGF_list = list(cgf1_alt1, cgf2_alt2))
#
#
#   K(c(0.1, 0.2, 0.3, 0.4, 0.5), c(12, 10, 0.64), concat_cgf) == K(c(0.1, 0.2, 0.3, 0.4, 0.5), c(12, 10, 0.64), my.cgf)
# }











