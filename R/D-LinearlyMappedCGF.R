# R/linearlyMappedCGF.R
# Object: linearlyMappedCGF






#' @keywords internal
.as_RTMB_mat <- function(A_) {
  # Already AD sparse
  if (inherits(A_, "adsparse")) return(A_)
  # Numeric sparse from Matrix
  if (inherits(A_, "sparseMatrix")) return(A_)
  # AD dense matrix: advector with a 'dim' attribute
  if (inherits(A_, "advector") && !is.null(attr(A_, "dim"))) return(A_)
  # Base numeric dense matrix
  if (is.matrix(A_)) return(Matrix::Matrix(A_, sparse = TRUE))
  stop("matrix_A(param) returned unsupported type: ", paste(class(A_), collapse = ", "),
       ". Expected numeric matrix, 'sparseMatrix', AD dense (advector with dim), or 'adsparse'.")
}







#' @keywords internal
.linearlyMappedCGF_internal <- function(cgf, matrix_A, ...){





  is_matrix_A_function <- is.function(matrix_A)
  A_fun <- NULL
  is_already_sparse <- FALSE

  if (!is_matrix_A_function) {
    # numeric (dense or sparse) at construction time --> standardize once
    if (!is.matrix(matrix_A) && !inherits(matrix_A, "sparseMatrix")) {
      stop("'matrix_A' must be a numeric matrix, a sparseMatrix, or a function returning one of these or an AD equivalent.")
    }
    if (!inherits(matrix_A, "sparseMatrix")) matrix_A <- Matrix::Matrix(matrix_A, sparse = TRUE)
    is_already_sparse <- TRUE
    A_fun <- function(param) matrix_A
  } else {
    ##### function case: DO NOT call Matrix::Matrix() on AD objects
    A_fun <- function(param) .as_RTMB_mat(matrix_A(param))
  }



  #### this is not needed, but I'll keep it for now; might be useful for exta checks (central spot for extra logic)
  get_sparse_A <- function(param) {
    if (is_already_sparse) matrix_A else A_fun(param)
  }

  child_K2_factor <- .K2_factor_method(cgf)
  extra_args <- list(...)
  extra_args <- extra_args[!vapply(extra_args, is.null, logical(1))]
  extra_names <- names(extra_args)
  unsupported_required <- intersect(
    extra_names,
    c("K", "K1", "K3operator", "K4operator")
  )
  if (length(unsupported_required) > 0L) {
    stop(
      "linearlyMappedCGF cannot override ",
      paste(unsupported_required, collapse = ", "),
      " through '...'; override the child CGF instead.",
      call. = FALSE
    )
  }
  contraction_pairs <- list(
    c("K4operatorAABB", "K4operatorAABB_factored"),
    c("K3K3operatorAABBCC", "K3K3operatorAABBCC_factored"),
    c("K3K3operatorABCABC", "K3K3operatorABCABC_factored")
  )
  factored_names <- vapply(contraction_pairs, `[[`, character(1), 2L)
  force_factored_contractions <- any(factored_names %in% extra_names)






  #---------------------------------------------
  # Single-block linearlyMapped CGF (iidReps = 1 OR {iidReps = NULL AND block_size = NULL})
  # Overrides for K, K1, K2, etc., where 'A_current' = get_sparse_A(parameter_vector)
  #---------------------------------------------

  # Key identity: K_Y(t) = K_X(A^T t) (with t assumed to be a column vector)
  K <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    cgf$K(as.vector(t(A_current) %*% tvec), parameter_vector)
  }

  # Key identity: K_Y' = A K_X'
  K1 <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    A_current %*% cgf$K1(as.vector(t(A_current) %*% tvec), parameter_vector)
  }

  # Key identity: K_Y'' = A K_X'' A^T
  K2 <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    k2_base <- cgf$K2(as.vector(t(A_current) %*% tvec), parameter_vector)
    A_current %*% k2_base %*% t(A_current)
  }

  # Key identity: K_Y(t) - t^T K_Y'(t) = K_X(A^T t) - t^T A K_X'(A^T t) = K_X(A^T t) - (A^T t)^T K_X'(A^T t)
  base_tilting_exponent <- cgf$.private_api$tilting_exponent
  tilting_exponent <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    base_tilting_exponent(as.vector(t(A_current) %*% tvec), parameter_vector)
  }

  # neg_ll: cgf's default neg_ll will be used (no override here).

  # Key identity: x^T K_Y'' y = x^T A K_X'' A^T y = (A^T x)^T K_X'' A^T y
  K2operator <- function(tvec, parameter_vector, x, y) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K2operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% x),
                   as.vector(t(A_current) %*% y)
    )
  }

  # Returns B K_Y'' B^T as a function of the supplied (non-parameter) argument B
  # Key identity: B K_Y'' B^T = B A K_X'' A^T B^T = (B A) K_X'' (B A)^T
  K2operatorAK2AT <- function(tvec, parameter_vector, B) {
    A_current <- get_sparse_A(parameter_vector)
    B_A <- B %*% A_current
    if (inherits(B_A, "denseMatrix") && !inherits(B_A, "adsparse")) {
      B_A <- as.matrix(B_A)
    }
    cgf$K2operatorAK2AT(as.vector(t(A_current) %*% tvec), parameter_vector, B_A)
  }

  K2_factor <- NULL
  if (!is.null(child_K2_factor)) {
    K2_factor <- function(tvec, parameter_vector, B) {
      A_current <- get_sparse_A(parameter_vector)
      B_A <- B %*% A_current
      if (inherits(B_A, "denseMatrix") && !inherits(B_A, "adsparse")) {
        B_A <- as.matrix(B_A)
      }
      child_K2_factor(
        as.vector(t(A_current) %*% tvec),
        parameter_vector,
        B_A
      )
    }
  }

  K3operator <- function(tvec, parameter_vector, v1, v2, v3) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K3operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% v1),
                   as.vector(t(A_current) %*% v2),
                   as.vector(t(A_current) %*% v3)
    )
  }

  K4operator <- function(tvec, parameter_vector, v1, v2, v3, v4) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K4operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% v1),
                   as.vector(t(A_current) %*% v2),
                   as.vector(t(A_current) %*% v3),
                   as.vector(t(A_current) %*% v4)
    )
  }


  # Factor Q before pulling it back.  The output-space Q used by the public
  # contractions is positive definite, whereas A^T Q A is only positive
  # semidefinite when A reduces dimension.  Passing A^T times a factor of Q to
  # the factored methods preserves the contraction without asking a child to
  # Cholesky-factor that singular pullback.
  factor_output_Q <- function(Q, normalize = TRUE) {
    if (inherits(Q, "Matrix") && !inherits(Q, "adsparse")) Q <- as.matrix(Q)
    chol_Q <- chol(Q)
    if (!normalize) {
      return(list(B = t(chol_Q), d = rep(1, nrow(chol_Q))))
    }
    diag_Q <- diag(chol_Q)
    list(
      B = t(chol_Q) %*% diag(
        1 / diag_Q,
        nrow = length(diag_Q),
        ncol = length(diag_Q)
      ),
      d = diag_Q * diag_Q
    )
  }

  mapped_Q_contraction <- function(self_object, tvec, parameter_vector, Q,
                                   dense_method, factored_name) {
    A_current <- get_sparse_A(parameter_vector)
    if (!force_factored_contractions && nrow(A_current) >= ncol(A_current)) {
      tA <- t(A_current)
      Q_inner <- tA %*% Q %*% A_current
      return(dense_method(
        as.vector(tA %*% tvec), parameter_vector, Q_inner
      ))
    }

    # Unit weights avoid an unnecessary dense column-normalization step here.
    # The factored contract only requires Q = B diag(d) B', which the raw
    # Cholesky factor satisfies exactly.
    Q_factor <- factor_output_Q(Q, normalize = FALSE)
    if (force_factored_contractions) {
      return(self_object$.private_api[[factored_name]](
        tvec, parameter_vector, Q_factor$B, Q_factor$d
      ))
    }

    tA <- t(A_current)
    cgf$.private_api[[factored_name]](
      as.vector(tA %*% tvec), parameter_vector,
      tA %*% Q_factor$B, Q_factor$d
    )
  }

  K4operatorAABB <- function(tvec, parameter_vector, Q) {
    mapped_Q_contraction(
      get("self", inherits = TRUE), tvec, parameter_vector, Q,
      cgf$K4operatorAABB, "K4operatorAABB_factored"
    )
  }

  K3K3operatorAABBCC <- function(tvec, parameter_vector, Q) {
    mapped_Q_contraction(
      get("self", inherits = TRUE), tvec, parameter_vector, Q,
      cgf$K3K3operatorAABBCC, "K3K3operatorAABBCC_factored"
    )
  }

  K3K3operatorABCABC <- function(tvec, parameter_vector, Q) {
    mapped_Q_contraction(
      get("self", inherits = TRUE), tvec, parameter_vector, Q,
      cgf$K3K3operatorABCABC, "K3K3operatorABCABC_factored"
    )
  }

  # Ordinary full-rank maps keep the established dense-Q contractions.  A
  # dimension-reducing map makes A' Q A singular; factor-capable compositions
  # also retain their stable thin representation.
  func_T <- function(tvec, parameter_vector) {
    self_object <- get("self", inherits = TRUE)
    private_api <- self_object$.private_api
    A_current <- get_sparse_A(parameter_vector)
    use_factored <- force_factored_contractions ||
      !is.null(.K2_factor_method(self_object)) ||
      nrow(A_current) < ncol(A_current)
    Q <- self_object$K2_solve(tvec, parameter_vector, diag(length(tvec)))
    if (inherits(Q, "Matrix") && !inherits(Q, "adsparse")) Q <- as.matrix(Q)

    if (!use_factored) {
      K3K3operatorABCABC_val <- self_object$K3K3operatorABCABC(
        tvec, parameter_vector, Q
      )
      K3K3operatorAABBCC_val <- self_object$K3K3operatorAABBCC(
        tvec, parameter_vector, Q
      )
      K4operatorAABB_val <- self_object$K4operatorAABB(
        tvec, parameter_vector, Q
      )
      return(
        K4operatorAABB_val / 8 -
          K3K3operatorAABBCC_val / 8 -
          K3K3operatorABCABC_val / 12
      )
    }

    Q_factor <- factor_output_Q(Q)
    B <- Q_factor$B
    d <- Q_factor$d

    K4_AABB <- private_api$K4operatorAABB_factored(tvec, parameter_vector, B, d)
    K3K3_AABBCC <- private_api$K3K3operatorAABBCC_factored(
      tvec, parameter_vector, B, d
    )
    K3K3_ABC <- private_api$K3K3operatorABCABC_factored(
      tvec, parameter_vector, B, d
    )
    K4_AABB / 8 - K3K3_AABBCC / 8 - K3K3_ABC / 12
  }

  # For the factored forms where Q = B D B^T and D has diagonal vector d, note that Q_inner = A^T Q A = (A^T B) D (A^T B)^T
  # Note about sizes: if A is n-by-m then B is n-by-r for some r, and A^T B is m-by-r
  base_K4operatorAABB_factored <- cgf$.private_api$K4operatorAABB_factored
  K4operatorAABB_factored <- function(tvec, parameter_vector, B, d) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B_inner <- tA %*% B
    base_K4operatorAABB_factored(as.vector(tA %*% tvec), parameter_vector, B_inner, d)
  }
  K4operatorAABB_factored <- .factored_delegate_mark(
    K4operatorAABB_factored,
    .factored_delegate_is_safe(base_K4operatorAABB_factored)
  )

  base_K3K3operatorAABBCC_factored <- cgf$.private_api$K3K3operatorAABBCC_factored
  K3K3operatorAABBCC_factored <- function(tvec, parameter_vector, B, d) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B_inner <- tA %*% B
    base_K3K3operatorAABBCC_factored(as.vector(tA %*% tvec), parameter_vector, B_inner, d)
  }

  base_K3K3operatorABCABC_factored <- cgf$.private_api$K3K3operatorABCABC_factored
  K3K3operatorABCABC_factored <- function(tvec, parameter_vector, B, d) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B_inner <- tA %*% B
    base_K3K3operatorABCABC_factored(as.vector(tA %*% tvec), parameter_vector, B_inner, d)
  }
  K3K3operatorAABBCC_factored <- .factored_delegate_mark(
    K3K3operatorAABBCC_factored,
    .factored_delegate_is_safe(base_K3K3operatorAABBCC_factored)
  )
  K3K3operatorABCABC_factored <- .factored_delegate_mark(
    K3K3operatorABCABC_factored,
    .factored_delegate_is_safe(base_K3K3operatorABCABC_factored)
  )

  # inequality constraints for the transformed variable Y = A * X are the same as those
  # for the original variable X, evaluated at the transformed input A.transpose() * tvec.
  ineq_constraint <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$ineq_constraint(as.vector(t(A_current) %*% tvec), parameter_vector)
  }


  rsim <- NULL
  if (isTRUE(cgf$has_rsim)) {
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      A_current <- get_sparse_A(parameter_vector)
      d_out <- nrow(A_current)
      d_in  <- ncol(A_current)

      # tilt propagation: t_x = A^T t_y
      t_inner <- if (is.null(tvec)) NULL else as.vector(t(A_current) %*% tvec)

      X <- cgf$rsim(
        n = n,
        vector_length = d_in,
        parameter_vector = parameter_vector,
        tvec = t_inner,
        flatten = FALSE,
        ...
      )

      as.matrix(A_current %*% X)  # d_out x n
    }
  }


  # # ------------------------------------------------------------------
  # # Optional simulator: if X ~ cgf and Y = A X then Y_sim = A %*% X_sim
  # # ------------------------------------------------------------------
  # simulate_fun <- NULL
  # if (isTRUE(cgf$has_rsim)) {
  #   simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
  #     A_current <- get_sparse_A(parameter_vector)
  #     d_out <- nrow(A_current)
  #     d_in <- ncol(A_current)
  #
  #     if (vector_length %% d_out != 0L) {
  #       stop("linearlyMappedCGF$rsim: 'vector_length' must be a multiple of nrow(A).", call. = FALSE)
  #     }
  #     blocks <- as.integer(vector_length %/% d_out)
  #
  #     out <- matrix(0, nrow = vector_length, ncol = n)
  #
  #     if (is.null(tvec)) {
  #       X <- cgf$rsim(
  #         n = n * blocks,
  #         vector_length = d_in,
  #         parameter_vector = parameter_vector,
  #         tvec = NULL,
  #         flatten = FALSE,
  #         ...
  #       )
  #       Y <- as.matrix(A_current %*% X)
  #       for (j in seq_len(n)) {
  #         cols <- ((j - 1L) * blocks + 1L):(j * blocks)
  #         out[, j] <- as.vector(Y[, cols, drop = FALSE])
  #       }
  #       return(out)
  #     }
  #
  #     tA <- t(A_current)
  #     for (b in seq_len(blocks)) {
  #       rows <- ((b - 1L) * d_out + 1L):(b * d_out)
  #       t_y  <- tvec[rows]
  #       t_x  <- as.vector(tA %*% t_y)
  #
  #       Xb <- cgf$rsim(
  #         n = n,
  #         vector_length = d_in,
  #         parameter_vector = parameter_vector,
  #         tvec = t_x,
  #         flatten = FALSE,
  #         ...
  #       )
  #       Yb <- as.matrix(A_current %*% Xb)
  #       out[rows, ] <- Yb
  #     }
  #
  #     out
  #   }
  # }

  # Build args list (names match createCGF parameters exactly)
  cgf_args <- list(
    K = K,
    K1 = K1,
    K2 = K2,
    K3operator = K3operator,
    K4operator = K4operator,
    ineq_constraint = ineq_constraint,
    analytic_tvec_hat = NULL,
    tilting_exponent = tilting_exponent,
    rsim = rsim,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    K2operator = K2operator,
    K2operatorAK2AT = K2operatorAK2AT,
    K2_factor = K2_factor,
    K2_factor_terminal = if (!is.null(K2_factor)) function() TRUE else NULL,
    op_name = c(cgf$call_history, "linearlyMappedCGF")
  )

  cgf_args$func_T <- func_T

  if (any(c("K2", "K2operatorAK2AT") %in% extra_names) &&
      !("K2_factor" %in% extra_names)) {
    cgf_args$K2_factor <- NULL
    cgf_args$K2_factor_terminal <- NULL
  }

  if ("K2" %in% extra_names) {
    for (dependent_name in c("K2operator", "K2operatorAK2AT")) {
      if (!(dependent_name %in% extra_names)) {
        cgf_args[[dependent_name]] <- NULL
      }
    }
  }

  for (pair in contraction_pairs) {
    if (pair[[1L]] %in% extra_names && !(pair[[2L]] %in% extra_names)) {
      cgf_args[[pair[[2L]]]] <- NULL
    }
  }

  do.call(createCGF, modifyList(cgf_args, extra_args))
}













#' @title CGF Object of a linearly mapped random variable \eqn{Y = A \, X}
#'
#' @description
#' Creates a CGF object for the random vector \eqn{Y = A(\theta) \, X}, where
#' \eqn{X} is described by the input CGF `cgf`. The argument `matrix_A` can be:
#'
#' - **A numeric matrix** (dense or sparse).
#' - **A function**: \eqn{\theta \mapsto A(\theta)}
#' See details for the possible options of `matrix_A`.
#'
#' If `matrix_A` is a function, it is called for each invocation of the CGF
#' methods to retrieve the current matrix (allowing parameter-dependent transformations).
#'
#' @details
#' Accepted types for \code{matrix_A}:
#' \itemize{
#'   \item \strong{numeric constant} matrix, dense or \code{Matrix} sparse. If dense, it is converted once
#'         at construction to a sparse \code{Matrix}.
#'   \item \strong{function} \eqn{\theta \mapsto A(\theta)} returning one of:
#'         \itemize{
#'           \item a numeric \code{matrix} (dense); it will be internally handled as sparse;
#'           \item a \code{Matrix} sparse matrix (preferred for large problems);
#'           \item an RTMB::AD-dense matrix (RTMB \code{advector} with a \code{dim} attribute);
#'           \item an RTMB::AD-sparse matrix (RTMB \code{adsparse}), created e.g. via
#'                 \code{A = RTMB::AD(Matrix::sparseMatrix(...)); A@x[] = ...}.
#'         }
#' }
#'
#'
#'
#'
#' @param cgf An object of class `CGF` for the base distribution \eqn{X}.
#' @param matrix_A Either a numeric matrix (dense or \code{Matrix} sparse), or a function
#'   \code{function(theta) -> A(theta)} returning one of: numeric dense matrix,
#'   \code{Matrix} sparse matrix, RTMB AD‑dense (an \code{advector} with a \code{dim})
#'   or RTMB \code{adsparse}.
#'
#' @param iidReps Either \code{"any"} (default) or a positive integer. See
#'   \code{\link{iidReplicatesCGF}} for the replication semantics.
#' @param ... Additional named arguments passed to `CGF` creation functions.
#'
#'
#' @examples
#' ## Example 1: constant numeric A
#' \dontrun{
#' lambda_fun <- function(theta) c(theta[1], theta[2])   # two Poisson rates
#' pois2 <- PoissonModelCGF(lambda = lambda_fun, iidReps = "any")
#'
#' A_const <- rbind(c(1, 0),
#'                  c(0.5, 1))
#' mapped  <- linearlyMappedCGF(cgf = pois2, matrix_A = A_const, iidReps = "any")
#'
#' B <- 3L # 3 replicated 2-d blocks
#' t_one <- c(0.10, -0.05)
#' y_one <- c(3.0,   4.0)
#' tvec  <- rep(t_one, B)
#' y <- rep(y_one, B)
#' theta <- c(2.0, 3.0)
#'
#' ## Identity check: K1_Y(t) = A K1_X(A^T t)
#' t_block <- t_one
#' A_now   <- A_const
#' lhs_K1  <- mapped$K1(t_block, theta)
#' poisk1 <- pois2$K1(as.vector(t(A_now) %*% t_block), theta)
#' rhs_K1  <- A_now %*% poisk1
#'
#' ## SPA negative log-likelihood + gradient + Hessian
#' res <- compute.spa.negll(
#'   parameter_vector = theta,
#'   observed.data    = y,
#'   cgf              = mapped,
#'   gradient         = TRUE,
#'   hessian          = TRUE,
#'   tvec_source      = "newton",
#'   spa_method       = "standard"
#' )
#' }
#' ## Example 2: A(theta) dense-vs-sparse (adsparse)
#' \dontrun{
#' library(Matrix)
#' library(RTMB)
#'
#' lambda_fun <- function(theta) c(theta[1], theta[2])   # base 2-d Poisson
#' pois2 <- PoissonModelCGF(lambda = lambda_fun, iidReps = "any")
#'
#' ## Dense A(theta) depending on theta[1]
#' A_theta_dense <- function(theta) {
#'   matrix(c(1, 0,
#'            0.5*theta[1], 1), 2, 2, byrow = TRUE)
#' }
#'
#' ## AD-sparse A(theta) with fixed sparsity pattern, AD-filled values
#' A_theta_sparse <- function(theta) {
#'   A0 <- sparseMatrix(i = c(1L, 2L, 2L),
#'                      j = c(1L, 1L, 2L),
#'                      x = c(1, 0, 1),
#'                      dims = c(2L, 2L))
#'   A  <- RTMB::AD(A0)
#'   A@x[] = c(1, 0.5*theta[1], 1)
#'   A
#' }
#'
#' mapped_dense  <- linearlyMappedCGF(cgf = pois2, matrix_A = A_theta_dense,  iidReps = "any")
#' mapped_sparse <- linearlyMappedCGF(cgf = pois2, matrix_A = A_theta_sparse, iidReps = "any")
#'
#' ## Data: B = 3 identical 2-vectors
#' B     <- 3L
#' y_one <- c(3.0, 4.0)
#' y     <- rep(y_one, B)
#' theta <- c(2.0, 3.0) # only theta[1] affects A(theta)
#'
#' nll_dense <- compute.spa.negll(theta, y, mapped_dense,
#'                                gradient=TRUE, hessian=TRUE,
#'                                tvec_source="newton", spa_method="standard")
#' nll_sparse <- compute.spa.negll(theta, y, mapped_sparse,
#'                                 gradient=TRUE, hessian=TRUE,
#'                                 tvec_source="newton", spa_method="standard")
#' }
#'
#' @return A `CGF` object for \eqn{Y = A \, X}.
#' @export
linearlyMappedCGF <- function(cgf, matrix_A, iidReps = "any", ...) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be an object inheriting from class 'CGF'.")
  .check_iidReps(iidReps)

  mapped <- .linearlyMappedCGF_internal(cgf, matrix_A, ...)


  block_size <- if (is.function(matrix_A)) {
    bs_fun <- function(param) {
      A_ <- matrix_A(param)
      as.integer(nrow(A_))
    }
    # only for printing
    attr(bs_fun, "label") <- "nrow(A(theta))" # cosmetic: used only for informative printing
    bs_fun
  } else {
    as.integer(nrow(matrix_A))
  }


  iidReplicatesCGF(cgf = mapped, iidReps = iidReps, block_size = block_size)
}
