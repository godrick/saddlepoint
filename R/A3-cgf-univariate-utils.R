# R/A2-cgf-univariate-utils.R
# Univariate CGF utilities:
# - "iidReps" is *always* either "any" or a positive integer (no NULL, no tolower).
# - Broadcasting is done with matrices (no lists) to keep AD happy.

# Validate iidReps in { "any", positive integer }
#' @noRd
.check_iidReps <- function(iidReps) {
  if (identical(iidReps, "any")) return(invisible())
  if (!is.numeric(iidReps) || length(iidReps) != 1 || !is.finite(iidReps) ||
      iidReps < 1 || iidReps != as.integer(iidReps)) {
    stop("'iidReps' must be \"any\" or a positive integer.")
  }
  invisible()
}

# Broadcast parameter *matrix* (M x D) to match length(tvec) under iidReps policy.
# Returns a matrix with N rows (N = length(tvec)) and D columns.
#' @noRd
.broadcast_univariate_matrix <- function(tvec, par_mat, iidReps) {
  if (!is.matrix(par_mat)) stop("Internal error: 'par_mat' must be a matrix.")
  M <- nrow(par_mat)
  if (M < 1L) stop("Parameter vectors must be non-empty.")
  N <- length(tvec)

  if (identical(iidReps, "any")) {
    if (N %% M != 0L) {
      stop(sprintf("'tvec' length %d is not a multiple of parameter length %d (iidReps=\"any\").", N, M))
    }
    reps <- N %/% M
  } else {
    reps <- as.integer(iidReps)
    expected <- M * reps
    if (N != expected) {
      stop(sprintf("'tvec' length %d != %d (= %d parameters * iidReps=%d).", N, expected, M, reps))
    }
  }

  if (reps == 1L && N == M) {
    par_mat
  } else {
    par_mat[rep.int(seq_len(M), times = reps), , drop = FALSE]
  }
}

# Build a univariate CGF object from *elementwise* formulas.
# The elementwise functions receive:
#   - tvec            : numeric/AD vector of length N
#   - par_mat_expanded: N x D matrix (already broadcast)
#
# split_param_to_mat(param) must return an M x D matrix of base parameters.
#' @noRd
.make_univariate_model_cgf_matrix <- function(
  K_elem, K1_elem, K2_elem, K3_elem, K4_elem, t_hat_elem,
  split_param_to_mat,
  iidReps,
  op_name,
  ineq_elem = NULL,
  rsim_elem = NULL,
  ...
) {
  .check_iidReps(iidReps)

  stopifnot(is.function(K_elem),
            is.function(K1_elem),
            is.function(K2_elem),
            is.function(K3_elem),
            is.function(K4_elem),
            is.function(t_hat_elem),
            is.function(split_param_to_mat))
  if (!is.null(ineq_elem)) stopifnot(is.function(ineq_elem))
  if (!is.null(rsim_elem)) stopifnot(is.function(rsim_elem))

  # Helper that prepares the expanded parameter matrix for a given vector 'vec'
  prep_par <- function(vec, param) {
    base_mat <- unname(split_param_to_mat(param))
    .broadcast_univariate_matrix(vec, base_mat, iidReps)
  }

  rsim_wrapper <- NULL
  if (!is.null(rsim_elem)) {
    rsim_wrapper <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      t_use <- if (is.null(tvec)) numeric(vector_length) else tvec
      pm <- prep_par(t_use, parameter_vector)
      rsim_elem(n = n, tvec = t_use, pm = pm, ...)
    }
  }

  createCGF_fromVectorisedFunctions(
    K_vectorized_func  = function(tvec, param) {
      pm <- prep_par(tvec, param)
      K_elem(tvec, pm)
    },
    K1_vectorized_func = function(tvec, param) {
      # print(class(tvec))
      # print(class(param))
      pm <- prep_par(tvec, param)
      # print(class(pm))
      K1_elem(tvec, pm)
    },
    K2_vectorized_func = function(tvec, param) {
      pm <- prep_par(tvec, param)
      K2_elem(tvec, pm)  # returns a length-N vector (interpreted as diag)
    },
    K3_vectorized_func = function(tvec, param) {
      pm <- prep_par(tvec, param)
      K3_elem(tvec, pm)
    },
    K4_vectorized_func = function(tvec, param) {
      pm <- prep_par(tvec, param)
      K4_elem(tvec, pm)
    },
    analytic_tvec_hat = function(x, param) {
      pm <- prep_par(x, param)
      t_hat_elem(x, pm)
    },
    ineq_constraint = if (!is.null(ineq_elem)) {
      function(tvec, param) {
        pm <- prep_par(tvec, param)
        ineq_elem(tvec, pm)
      }
    } else NULL,
    rsim = rsim_wrapper,
    op_name = op_name,
    ...
  )
}
