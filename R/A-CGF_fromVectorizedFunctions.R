# R/A-CGF_fromVectorizedFunctions.R
# Objects: VectorizedFunctionsCGF, createCGFfromVectorizedFunctions, createCGF_fromVectorisedFunctions

#' A Vectorized CGF Class
#'
#' @description
#' `VectorizedFunctionsCGF` inherits from `CGF` and makes it easier to create a CGF
#' when you have vectorized functions.
#'
#' This class constructs the necessary `K()`, `K1()`, `K2()`, `K3operator()`,
#' and `K4operator()` from your vectorized forms. Optional methods are resolved
#' by starting from vectorized defaults and then merging user-supplied overrides.
#'
#' @noRd
VectorizedFunctionsCGF <- R6::R6Class(
  "VectorizedFunctionsCGF",
  inherit = CGF,

  private = list(
    K_vectorized  = NULL,
    K1_vectorized = NULL,
    K2_vectorized = NULL,
    K3_vectorized = NULL,
    K4_vectorized = NULL
  ),

  public = list(
    initialize = function(
      K_vectorized,
      K1_vectorized,
      K2_vectorized,
      K3_vectorized,
      K4_vectorized,
      ineq_constraint = NULL,
      analytic_tvec_hat = NULL,
      rsim = NULL,
      op_name = "UnnamedOperation",
      tilting_exponent = NULL,
      neg_ll = NULL,
      func_T = NULL,
      K4operatorAABB = NULL,
      K3K3operatorAABBCC = NULL,
      K3K3operatorABCABC = NULL,
      K4operatorAABB_factored = NULL,
      K3K3operatorAABBCC_factored = NULL,
      K3K3operatorABCABC_factored = NULL,
      K2operator = NULL,
      K2operatorAK2AT = NULL,
      ...
    ) {
      if (!is.function(K_vectorized) || !is.function(K1_vectorized) ||
          !is.function(K2_vectorized) || !is.function(K3_vectorized) ||
          !is.function(K4_vectorized)) {
        stop(
          "K_vectorized, K1_vectorized, K2_vectorized, K3_vectorized, and K4_vectorized must all be functions.",
          call. = FALSE
        )
      }

      private$K_vectorized  <- K_vectorized
      private$K1_vectorized <- K1_vectorized
      private$K2_vectorized <- K2_vectorized
      private$K3_vectorized <- K3_vectorized
      private$K4_vectorized <- K4_vectorized

      # fallback implementations based on the vectorized K/K1/K2/K3/K4
      default_methods <- list(
        tilting_exponent = function(tvec, p) {
          K_vals <- private$K_vectorized(tvec, p)
          K1_vals <- private$K1_vectorized(tvec, p)
          sum(K_vals - tvec * K1_vals)
        },
        neg_ll = function(tvec, p) {
          K2_vals <- private$K2_vectorized(tvec, p)
          K_vals <- private$K_vectorized(tvec, p)
          K1_vals <- private$K1_vectorized(tvec, p)
          tilting_vals <- K_vals - tvec * K1_vals
          sum(0.5 * log(2 * pi * K2_vals) - tilting_vals)
        },
        K2operatorAK2AT = function(tvec, p, A) {
          K2_vals <- private$K2_vectorized(tvec, p)
          A %*% (K2_vals * t(A))
        },
        func_T = function(tvec, p) {
          k2val <- private$K2_vectorized(tvec, p)
          k2sq_val <- k2val * k2val
          k3val <- private$K3_vectorized(tvec, p)
          k4val <- private$K4_vectorized(tvec, p)
          sum(k4val/(8 * k2sq_val) - 5*(k3val*k3val)/(24 * k2sq_val * k2val))
        },
        K4operatorAABB = function(tvec, p, Q) {
          sum(private$K4_vectorized(tvec, p) * diag(Q) * diag(Q))
        },
        K3K3operatorAABBCC = function(tvec, p, Q) {
          k3_vals <- private$K3_vectorized(tvec, p)
          sum((diag(Q) * k3_vals) %*% Q %*% (diag(Q) * k3_vals))
        },
        K3K3operatorABCABC = function(tvec, p, Q) {
          k3_vals <- private$K3_vectorized(tvec, p)
          mat_k3_vals <- diag(k3_vals, nrow = length(tvec))
          sum(mat_k3_vals %*% (Q * Q * Q) %*% mat_k3_vals)
        },
        K4operatorAABB_factored = function(tvec, p, A, d) {
          diag_Q <- as.vector(rowSums(A * t(d * t(A))))
          sum(private$K4_vectorized(tvec, p) * diag_Q * diag_Q)
        },
        K3K3operatorAABBCC_factored = function(tvec, p, A, d) {
          k3_vals <- private$K3_vectorized(tvec, p)
          balanced <- .balance_factored_Q(
            A, d, "vectorized K3K3operatorAABBCC_factored"
          )
          A <- balanced$A
          d <- balanced$d
          diag_Q <- as.vector(rowSums(A * t(d * t(A))))
          z <- as.vector(crossprod(A, k3_vals * diag_Q))
          sum((d * z) * z)
        },
        K3K3operatorABCABC_factored = function(tvec, p, A, d) {
          k3_vals <- private$K3_vectorized(tvec, p)
          r <- length(d)
          if (r == 0L) return(0 * sum(k3_vals))

          # Factor-space contraction is O(n r^3) and avoids an n-by-n Q when
          # the rank is genuinely small.  A dense elementwise contraction is
          # cheaper when r is large, so select by the leading operation counts.
          if (as.double(r) * r <= nrow(A)) {
            balanced <- .balance_factored_Q(
              A, d, "vectorized K3K3operatorABCABC_factored"
            )
            A <- balanced$A
            d <- balanced$d
            d_outer <- tcrossprod(d)
            total <- 0 * sum(k3_vals)
            for (p_index in seq_len(r)) {
              C_slice <- crossprod(
                A,
                A * as.vector(k3_vals * A[, p_index])
              )
              total <- total +
                d[p_index] * sum(d_outer * (C_slice * C_slice))
            }
            return(total)
          }

          Q <- A %*% (d * t(A))
          sum(tcrossprod(k3_vals) * (Q * Q * Q))
        }
      )

      ## here we collect named optional args from the function signature
      optional_overrides <- list(
        tilting_exponent = tilting_exponent,
        neg_ll = neg_ll,
        func_T = func_T,
        K4operatorAABB = K4operatorAABB,
        K3K3operatorAABBCC = K3K3operatorAABBCC,
        K3K3operatorABCABC = K3K3operatorABCABC,
        ineq_constraint = ineq_constraint,
        analytic_tvec_hat = analytic_tvec_hat,
        rsim = rsim,
        K4operatorAABB_factored = K4operatorAABB_factored,
        K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
        K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
        K2operator = K2operator,
        K2operatorAK2AT = K2operatorAK2AT
      )
      optional_overrides <- optional_overrides[
        !vapply(optional_overrides, is.null, logical(1))
      ]

      ## collect additional named overrides; it filters NULLs too and warns on unnamed entries
      extra_methods <- list(...)
      if (length(extra_methods) > 0L && is.null(names(extra_methods))) {
        warning("Unnamed entries in '...' are ignored. Please provide named overrides.", call. = FALSE)
        extra_methods <- list()
      }
      if (length(extra_methods) > 0L) {
        extra_methods <- extra_methods[
          !vapply(extra_methods, is.null, logical(1))
        ]
      }

      ## merge in precedence order; default_methods (lowest) --> explicit args like neg_ll=... --> ... (highest)
      user_methods <- modifyList(optional_overrides, extra_methods)
      user_names <- names(user_methods)
      correction_methods <- c(
        "K4operatorAABB",
        "K3K3operatorAABBCC",
        "K3K3operatorABCABC",
        "K4operatorAABB_factored",
        "K3K3operatorAABBCC_factored",
        "K3K3operatorABCABC_factored"
      )
      if (!("func_T" %in% user_names) &&
          any(correction_methods %in% user_names)) {
        default_methods$func_T <- NULL
      }
      if (!("neg_ll" %in% user_names) &&
          any(c("logdetK2", "tilting_exponent") %in% user_names)) {
        default_methods$neg_ll <- NULL
      }
      # Only explicit dense overrides need a compatibility bridge.  Bridging
      # the built-in dense defaults would reconstruct a large dense Q and erase
      # the thin factor supplied by a dimension-reducing composition.
      user_methods <- .add_factored_contraction_bridges(user_methods)
      resolved_methods <- modifyList(default_methods, user_methods)
      ## resolved_methods is passed into super$initialize

      ### effect: explicit args stay user-facing, but (...) still gives extensibility and can //??intentionally// override anything


      init_args <- c(
        list(
          K = function(tvec, p) { sum(private$K_vectorized(tvec, p)) },
          K1 = function(tvec, p) { private$K1_vectorized(tvec, p) },
          K2 = function(tvec, p) { diag(private$K2_vectorized(tvec, p), nrow = length(tvec)) },
          K3operator = function(tvec, p, v1, v2, v3) { sum(private$K3_vectorized(tvec, p) * v1 * v2 * v3) },
          K4operator = function(tvec, p, v1, v2, v3, v4) { sum(private$K4_vectorized(tvec, p) * v1 * v2 * v3 * v4) },
          op_name = op_name
        ),
        resolved_methods
      )

      do.call(super$initialize, init_args)
    }
  )
)


#' Create a `CGF` object from vectorized functions
#'
#' @description
#' This function allows you to create a `CGF` object using vectorized functions
#' along with any optional operators or methods.
#'
#' @param K_vectorized A function of the form \code{function(tvec, param) -> numeric vector} that returns the CGF values.
#' @param K1_vectorized A function of the form \code{function(tvec, param) -> numeric vector} that returns the first derivative values.
#' @param K2_vectorized A function of the form \code{function(tvec, param) -> numeric vector} that returns the second derivative values.
#' @param K3_vectorized,K4_vectorized Similar vectorized functions for the third and fourth derivatives.
#' @param ineq_constraint Optional inequality constraint function.
#' @param analytic_tvec_hat Optional `tvec_hat` function override.
#' @param op_name Optional character string indicating the name of the operation or transformation being performed.
#' @param tilting_exponent Optional tilting exponent function override.
#' @param neg_ll Optional neg_ll function override.
#' @param func_T Optional func_T function override.
#' @param rsim Optional simulation method. A function of the form
#'   \code{function(n, vector_length, parameter_vector, tvec = NULL, ...)} returning
#'   a numeric vector or matrix of length \code{n * vector_length} or a
#'   \code{vector_length x n} matrix.
#' @param K4operatorAABB,K3K3operatorAABBCC,K3K3operatorABCABC Optional operator overrides.
#' @param K4operatorAABB_factored,K3K3operatorAABBCC_factored,K3K3operatorABCABC_factored Optional factored operator overrides.
#' @param K2operator,K2operatorAK2AT Optional operator overrides.
#' @param ... Any additional named optional methods.
#'
#' @return A `CGF` object.
#' @export
createCGFfromVectorizedFunctions <- function(
  K_vectorized,
  K1_vectorized,
  K2_vectorized,
  K3_vectorized,
  K4_vectorized,
  ineq_constraint = NULL,
  analytic_tvec_hat = NULL,
  op_name = "UnnamedOperation",
  tilting_exponent = NULL,
  neg_ll = NULL,
  func_T = NULL,
  rsim = NULL,
  K4operatorAABB = NULL,
  K3K3operatorAABBCC = NULL,
  K3K3operatorABCABC = NULL,
  K4operatorAABB_factored = NULL,
  K3K3operatorAABBCC_factored = NULL,
  K3K3operatorABCABC_factored = NULL,
  K2operator = NULL,
  K2operatorAK2AT = NULL,
  ...
) {
  do.call(VectorizedFunctionsCGF$new, c(
    list(
      K_vectorized = K_vectorized,
      K1_vectorized = K1_vectorized,
      K2_vectorized = K2_vectorized,
      K3_vectorized = K3_vectorized,
      K4_vectorized = K4_vectorized,
      ineq_constraint = ineq_constraint,
      analytic_tvec_hat = analytic_tvec_hat,
      op_name = op_name,
      tilting_exponent = tilting_exponent,
      neg_ll = neg_ll,
      func_T = func_T,
      rsim = rsim,
      K4operatorAABB = K4operatorAABB,
      K3K3operatorAABBCC = K3K3operatorAABBCC,
      K3K3operatorABCABC = K3K3operatorABCABC,
      K4operatorAABB_factored = K4operatorAABB_factored,
      K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
      K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
      K2operator = K2operator,
      K2operatorAK2AT = K2operatorAK2AT
    ),
    list(...)
  ))
}

#' @rdname createCGFfromVectorizedFunctions
#' @export
createCGF_fromVectorisedFunctions <- function(...) {
  createCGFfromVectorizedFunctions(...)
}
