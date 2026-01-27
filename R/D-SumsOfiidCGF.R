# ------------------------------------------------------------------
# R/D-SumsOfiidCGF.R
# Main function: sumOfiidCGF
#
# Purpose
# -------
# Create a CGF object for the sum of n i.i.d. copies of a base CGF:
#
#   Y = X_1 + ... + X_n,
#
# where the X_i are i.i.d. with CGF K_X(t;theta).
#
# Then:
#   K_Y(t;theta) = n * K_X(t;theta).
#
# Replication semantics
# ---------------------
# n is the number of summands INSIDE each Y.
# The arguments `iidReps` and `block_size` (if supplied) refer to i.i.d.
# replication of the *resulting variable Y* across observations.


#' @keywords internal
.scale_like_K2 <- function(K2_val, scalar_) {
  # Scale a K2-like object by a scalar while preserving matrix-ness for
  # RTMB AD-dense matrices (advector with dim attribute).

  if (length(K2_val) == 1) {
    # Scalar Hessian: force 1x1 matrix (helps downstream code that expects a matrix).
    return(scalar_ %*% K2_val)
  }

  out <- K2_val * scalar_

  # RTMB AD-dense matrices may lose the 'dim' attribute under scalar multiplication.
  if (!is.null(dim(K2_val)) && is.null(dim(out))) {
    attr(out, "dim") <- dim(K2_val)
  }

  out
}


#' @keywords internal
.sumOfiidCGF_internal <- function(cgf, n_fn, ...) {

  # private methods from the base CGF
  base_tilt <- cgf$.private_api$tilting_exponent
  base_T    <- cgf$.private_api$func_T

  .get_n <- function(param) {
    n_val <- n_fn(param)
    if (length(n_val) != 1L) stop("sumOfiidCGF: n(theta) must return a scalar.")
    # Only enforce positivity when n is plain numeric
    if (!inherits(n_val, "advector") && (!is.finite(n_val) || n_val <= 0)) {
      stop("sumOfiidCGF: n(theta) must be a finite scalar > 0.")
    }
    n_val
  }



  Kfun  <- function(tvec, param) {
    n_val <- .get_n(param)
    n_val * cgf$K(tvec, param)
  }

  K1fun <- function(tvec, param) {
    n_val <- .get_n(param)
    n_val * cgf$K1(tvec, param)
  }

  K2fun <- function(tvec, param) {
    n_val <- .get_n(param)
    .scale_like_K2(cgf$K2(tvec, param), n_val)
  }

  K3opfun <- function(tvec, param, v1, v2, v3) {
    n_val <- .get_n(param)
    n_val * cgf$K3operator(tvec, param, v1, v2, v3)
  }

  K4opfun <- function(tvec, param, v1, v2, v3, v4) {
    n_val <- .get_n(param)
    n_val * cgf$K4operator(tvec, param, v1, v2, v3, v4)
  }


  # tilting_exponent_Y(t) = n * tilting_exponent_X(t)
  tiltingfun <- function(tvec, param) {
    n_val <- .get_n(param)
    n_val * base_tilt(tvec, param)
  }


  # ------------------------------------------------------------------
  # Hessian utilities (important for speed + avoids scalar-matrix pitfalls)
  # ------------------------------------------------------------------

  # K2_Y = n * K2_X  =>  (K2_Y)^{-1} rhs = (1/n) * (K2_X)^{-1} rhs
  K2_solve_fun <- function(tvec, param, rhs) {
    n_val <- .get_n(param)
    cgf$K2_solve(tvec, param, rhs) / n_val
  }

  # log det(K2_Y) = log det(n*K2_X) = log det(K2_X) + d * log(n)
  # where d = dim(tvec) for this evaluation.
  logdetK2_fun <- function(tvec, param) {
    n_val <- .get_n(param)
    cgf$logdetK2(tvec, param) + length(tvec) * log(n_val)
  }


  # ------------------------------------------------------------------
  # Bilinear/linear operators
  # ------------------------------------------------------------------

  K2opfun <- function(tvec, param, x, y) {
    n_val <- .get_n(param)
    n_val * cgf$K2operator(tvec, param, x, y)
  }

  K2opAK2ATfun <- function(tvec, param, A) {
    n_val <- .get_n(param)
    .scale_like_K2(cgf$K2operatorAK2AT(tvec, param, A), n_val)
  }


  # ------------------------------------------------------------------
  # Higher-order correction term
  # ------------------------------------------------------------------

  # For Y = sum_{i=1}^n X_i:
  #   K^{(r)}_Y = n K^{(r)}_X,  r>=1
  #   Q_Y = (K2_Y)^{-1} = (1/n) Q_X
  # and the standard SPA correction term T scales as 1/n:
  #   T_Y = T_X / n.
  func_Tfun <- function(tvec, param) {
    n_val <- .get_n(param)
    base_T(tvec, param) / n_val
  }


  # ------------------------------------------------------------------
  # Higher-order Q-operators
  # ------------------------------------------------------------------

  K4AABBfun <- function(tvec, param, Q1, Q2) {
    n_val <- .get_n(param)
    n_val * cgf$K4operatorAABB(tvec, param, Q1, Q2)
  }

  K3K3AABBCCfun <- function(tvec, param, Q1, Q2, Q3) {
    n_val <- .get_n(param)
    (n_val * n_val) * cgf$K3K3operatorAABBCC(tvec, param, Q1, Q2, Q3)
  }

  K3K3ABCABCfun <- function(tvec, param, Q1, Q2, Q3) {
    n_val <- .get_n(param)
    (n_val * n_val) * cgf$K3K3operatorABCABC(tvec, param, Q1, Q2, Q3)
  }


  # ------------------------------------------------------------------
  # Factored higher-order operators
  # ------------------------------------------------------------------

  # These reuse the BASE CGF's private factored methods (if they are optimized),
  # then apply the appropriate scaling.

  base_K4AABB_factored   <- cgf$.private_api$K4operatorAABB_factored
  base_K3K3AABBCC_fact   <- cgf$.private_api$K3K3operatorAABBCC_factored
  base_K3K3ABCABC_fact   <- cgf$.private_api$K3K3operatorABCABC_factored

  K4AABB_factored_fun <- function(tvec, param, A1, d1, A2, d2) {
    n_val <- .get_n(param)
    n_val * base_K4AABB_factored(tvec, param, A1, d1, A2, d2)
  }

  K3K3AABBCC_factored_fun <- function(tvec, param, A1, d1, A2, d2, A3, d3) {
    n_val <- .get_n(param)
    (n_val * n_val) * base_K3K3AABBCC_fact(tvec, param, A1, d1, A2, d2, A3, d3)
  }

  K3K3ABCABC_factored_fun <- function(tvec, param, A1, d1, A2, d2, A3, d3) {
    n_val <- .get_n(param)
    (n_val * n_val) * base_K3K3ABCABC_fact(tvec, param, A1, d1, A2, d2, A3, d3)
  }


  # ------------------------------------------------------------------
  # Constraints + analytic t-hat
  # ------------------------------------------------------------------

  ineqfun <- function(tvec, param) cgf$ineq_constraint(tvec, param)

  # Analytic t-hat mapping:
  #   Solve n*K1_X(t) = y  <=>  K1_X(t) = y/n.
  if (isTRUE(cgf$has_analytic_tvec_hat)) {
    analytic_tvec_hat_func <- function(y, param) {
      n_val <- .get_n(param)
      cgf$analytic_tvec_hat(y / n_val, param)
    }
  } else {
    analytic_tvec_hat_func <- NULL
  }


  # simulation (only if base cgf can simulate)
  simulate_fun <- NULL
  if (isTRUE(cgf$has_simulate)) {
    simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      n_num <- .get_n(parameter_vector)
      if (abs(n_num - round(n_num)) > 1e-8) {
        stop("sumOfiidCGF$rsim: n(theta) must be an integer to simulate a sum of i.i.d. terms.",
             call. = FALSE)
      }
      n_int <- as.integer(round(n_num))

      if (n_int == 1L) {
        return(cgf$rsim(
          n = n,
          vector_length = vector_length,
          parameter_vector = parameter_vector,
          tvec = tvec,
          flatten = FALSE,
          ...
        ))
      }

      X_all <- cgf$rsim(
        n = n * n_int,
        vector_length = vector_length,
        parameter_vector = parameter_vector,
        tvec = tvec,
        flatten = FALSE,
        ...
      )

      out <- matrix(0, nrow = vector_length, ncol = n)
      pos <- 1L
      for (j in seq_len(n)) {
        out[, j] <- rowSums(X_all[, pos:(pos + n_int - 1L), drop = FALSE])
        pos <- pos + n_int
      }
      out
    }
  }



  # ------------------------------------------------------------------
  # Build the new CGF
  # ------------------------------------------------------------------

  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,

    tilting_exponent = tiltingfun,
    func_T           = func_Tfun,

    K2_solve  = K2_solve_fun,
    logdetK2  = logdetK2_fun,
    rsim = simulate_fun,

    K2operator      = K2opfun,
    K2operatorAK2AT = K2opAK2ATfun,

    K4operatorAABB     = K4AABBfun,
    K3K3operatorAABBCC = K3K3AABBCCfun,
    K3K3operatorABCABC = K3K3ABCABCfun,

    K4operatorAABB_factored     = K4AABB_factored_fun,
    K3K3operatorAABBCC_factored = K3K3AABBCC_factored_fun,
    K3K3operatorABCABC_factored = K3K3ABCABC_factored_fun,

    ineq_constraint = ineqfun,
    analytic_tvec_hat = analytic_tvec_hat_func,

    op_name = c(cgf$call_history, "sumOfiidCGF"),
    ...
  )
}




#' Create a CGF object for the sum of `n` i.i.d. random variables.
#'
#' @description
#' Given a base CGF object `cgf` describing a random variable \eqn{X},
#' this function returns a new CGF object for the sum of `n` i.i.d. copies of \eqn{X},
#' i.e. \eqn{Y = X_1 + \cdots + X_n}.
#'
#' @details
#' \deqn{K_Y(t;\theta) = n(\theta)\,K_X(t;\theta).}
#'
#' For every order \eqn{r \ge 1},
#' \deqn{K_Y^{(r)}(t;\theta) = n(\theta)\,K_X^{(r)}(t;\theta).}
#'
#' \itemize{
#'   \item `n` is the number of summands inside each \eqn{Y}.
#'   \item `iidReps`/`block_size` (if supplied) describe i.i.d. replication of \eqn{Y}
#'         across observations \eqn{Y_1,\ldots,Y_B}.
#' }
#'
#' **Replication (iidReps / block_size):**
#' \itemize{
#'   \item If both `iidReps` and `block_size` are `NULL`, no replication is applied.
#'   \item If `block_size` is provided but `iidReps` is `NULL`, we set `iidReps = "any"` and
#'         infer the number of blocks from `length(tvec) / block_size` at evaluation time.
#'   \item If `iidReps = "any"`, then `block_size` must be provided.
#'   \item If `iidReps` is a positive integer, `block_size` may be `NULL`, though providing it is safer.
#' }
#'
#'
#' @param cgf A `CGF` object describing the base random variable \eqn{X}.
#' @param n Either a positive scalar, or a function/adaptor mapping \eqn{\theta \mapsto n(\theta)}.
#' @param block_size Optional. Positive integer giving the dimension of one observation of \eqn{Y}.
#' @param iidReps Optional. Either `NULL`, `"any"`, or a positive integer.
#' @param ... Additional arguments passed to CGF creation.
#'
#' @examples
#' \dontrun{
#' base_cgf <- NormalCGF
#' sum5 <- sumOfiidCGF(base_cgf, n=5)
#' theta <- c(1.2, 0.8)
#' tvec  <- c(0.1, -0.2, 0.0)
#'
#' all.equal(sum5$K(tvec, theta), 5*base_cgf$K(tvec, theta))
#' all.equal(sum5$K1(tvec, theta), 5*base_cgf$K1(tvec, theta))
#' }
#'
#' @return A `CGF` object for \eqn{Y = \sum_{i=1}^{n(\theta)} X_i}.
#' @export
sumOfiidCGF <- function(cgf,
                        n,
                        block_size = NULL,
                        iidReps = NULL,
                        ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be an object inheriting from class 'CGF'.")

  # n(theta)
  if (is.numeric(n)) {
    if (length(n) != 1L || is.na(n) || n <= 0) stop("'n' must be a positive numeric scalar.")
    n_const <- as.numeric(n)
    # Keep RTMB-friendly behavior: if theta is AD, n becomes AD-constant.
    n_fn <- function(theta) {
      if (length(theta) >= 1L) return(n_const + 0*theta[1])
      n_const
    }
  } else {
    n_fn <- validate_function_or_adaptor(n)
  }

  base_cgf <- .sumOfiidCGF_internal(cgf, n_fn, ...)

  # No explicit replication requested.
  if (is.null(block_size) && is.null(iidReps)) return(base_cgf)

  # Interpret block_size-only as iidReps = "any" (infer number of blocks later).
  if (is.null(iidReps)) iidReps <- "any"

  .check_iidReps(iidReps)

  # Explicit single replicate.
  if (is.numeric(iidReps) && iidReps == 1L && is.null(block_size)) return(base_cgf)

  if (identical(iidReps, "any") && is.null(block_size)) {
    stop("sumOfiidCGF(): iidReps='any' requires a non-NULL 'block_size'.")
  }

  iidReplicatesCGF(cgf = base_cgf, iidReps = iidReps, block_size = block_size)
}
