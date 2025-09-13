# R/C-PoissonCGF.R
# Objects: PoissonCGF, PoissonModelCGF

# --- Elementwise formulas: used by both ready CGF and model CGF -------------
# par_mat columns: [1] lambda
.poisson_K      <- function(tvec, pm) { pm[,1] * (exp(tvec) - 1) }
.poisson_K1     <- function(tvec, pm) { pm[,1] *  exp(tvec)       }
.poisson_K2     <- function(tvec, pm) { pm[,1] *  exp(tvec)       }
.poisson_K3     <- function(tvec, pm) { pm[,1] *  exp(tvec)       }
.poisson_K4     <- function(tvec, pm) { pm[,1] *  exp(tvec)       }
.poisson_That   <- function(x,    pm) { log(x / pm[,1]) }

# Split parameter vector into an Mx1 matrix (no checks yet; general on purpose)
.poisson_split_param_to_mat <- function(param) {
  # Accepts scalar or vector lambda; returns M x 1 matrix.
  matrix(param, ncol = 1L)
}

#' Poisson CGF Object
#'
#' Ready-to-use CGF for Poisson. Accepts scalar or vector `lambda`.
#' If `length(tvec)` is a multiple of `length(lambda)`, evaluation proceeds (iidReps="any").
#'
#' @examples
#' # Evaluate K at t = 0.1 for lambda = 2
#' # PoissonCGF$K(0.1, 2)
#'
#' @export
PoissonCGF <- .make_univariate_model_cgf_matrix(
  K_elem = .poisson_K, K1_elem = .poisson_K1, K2_elem = .poisson_K2,
  K3_elem = .poisson_K3, K4_elem = .poisson_K4, That_elem = .poisson_That,
  split_param_to_mat = .poisson_split_param_to_mat,
  iidReps = "any",                     # always "any" by default
  op_name = "PoissonCGF"
)

#' Create a Parametric Poisson CGF Object
#'
#' @description
#' Poisson with rate(s) `lambda(theta)`. If `iidReps = "any"` (default),
#' `length(tvec)` must be a multiple of `length(lambda(theta))`. If `iidReps = m`,
#' then `length(tvec)` must be `m * length(lambda(theta))`.
#'
#' @param lambda A function or `adaptor` mapping `theta` -> scalar or vector of rates.
#' @param iidReps Either `"any"` or a positive integer. Default `"any"`.
#' @param ... Passed through to the CGF creation (e.g., to override operators).
#' @return A `CGF` object.
#'
#' @examples
#' # ex 1: i.i.d. scenario for univariate Poisson(rate = 2)
#' ex_iid <- PoissonModelCGF(lambda = function(x) x[1])
#' # OR ex_iid <- PoissonModelCGF(lambda = function(x) x[1], iidReps = 3)
#' ex_iid$K1(rep(0,3), 2)
#'
#' # ex 2: non-identical scenario with replication: lambda returns c(2,5), iidReps=3
#' ex_repeat <- PoissonModelCGF(lambda = adaptor(indices = 1:2), iidReps = 3)
#' # OR ex_repeat <- PoissonModelCGF(lambda = function(x) c(2,5) )
#' ex_repeat$K1(rep(0,6), c(2,5))
#'
#' @export
PoissonModelCGF <- function(lambda, iidReps = "any", ...) {
  .check_iidReps(iidReps)  # only "any" or positive integer

  lambda_fn <- validate_function_or_adaptor(lambda)

  base <- .make_univariate_model_cgf_matrix(
    K_elem = .poisson_K, K1_elem = .poisson_K1, K2_elem = .poisson_K2,
    K3_elem = .poisson_K3, K4_elem = .poisson_K4, That_elem = .poisson_That,
    split_param_to_mat = .poisson_split_param_to_mat,
    iidReps = iidReps,                   # enforce "any" or a fixed replication
    op_name = "PoissonModelCGF",
    ...
  )

  adaptCGF(cgf = base, adaptor = lambda_fn)
}
