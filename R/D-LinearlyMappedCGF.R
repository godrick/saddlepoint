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






  #---------------------------------------------
  # # Single-block linearlyMapped CGF (iidReps = 1 OR {iidReps = NULL AND block_size = NULL})
  # # Overrides for K, K1, K2, etc., where 'A_current' = get_sparse_A(parameter_vector)
  #---------------------------------------------

  # Key identity: K_Y(t) = K_X(A^T t) (with t assumed to be a column vector)
  Kfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    cgf$K(t(A_current) %*% tvec, parameter_vector)
  }

  # Key identity: K_Y' = A K_X'
  K1fun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    if (nrow(A_current) != length(tvec)) stop("Dimension mismatch: nrow(matrix_A) != length(tvec).")
    A_current %*% cgf$K1(as.vector(t(A_current) %*% tvec), parameter_vector)
  }

  # Key identity: K_Y'' = A K_X'' A^T
  K2fun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    k2_base <- cgf$K2(as.vector(t(A_current) %*% tvec), parameter_vector)
    A_current %*% k2_base %*% t(A_current)
  }

  # Key identity: K_Y(t) - t^T K_Y'(t) = K_X(A^T t) - t^T A K_X'(A^T t) = K_X(A^T t) - (A^T t)^T K_X'(A^T t)
  tilting_exponent <- cgf$.get_private_method("tilting_exponent")
  tiltingfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    tilting_exponent(as.vector(t(A_current) %*% tvec), parameter_vector)
  }

  # neg_ll: cgf's neg_ll will be used.
  # # neg_ll <- cgf$.get_private_method("neg_ll")
  # negllfun <- NULL

  # Key identity: x^T K_Y'' y = x^T A K_X'' A^T y = (A^T x)^T K_X'' A^T y
  K2operatorfun <- function(tvec, parameter_vector, x, y) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K2operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% x),
                   as.vector(t(A_current) %*% y),
    )
  }

  # Returns B K_Y'' B^T as a function of the supplied (non-parameter) argument B
  # Key identity: B K_Y'' B^T = B A K_X'' A^T B^T = (B A) K_X'' (B A)^T
  K2operatorAK2ATfun <- function(tvec, parameter_vector, B) {
    A_current <- get_sparse_A(parameter_vector)
    B_A <- B %*% A_current
    # cgf$K2operatorAK2AT(as.vector(t(A_current) %*% tvec), parameter_vector, B_A) %*% t(B_A)
    cgf$K2operatorAK2AT(as.vector(t(A_current) %*% tvec), parameter_vector, B_A)
  }

  K3operatorfun <- function(tvec, parameter_vector, v1, v2, v3) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K3operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% v1),
                   as.vector(t(A_current) %*% v2),
                   as.vector(t(A_current) %*% v3)
    )
  }

  K4operatorfun <- function(tvec, parameter_vector, v1, v2, v3, v4) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$K4operator(as.vector(t(A_current) %*% tvec),
                   parameter_vector,
                   as.vector(t(A_current) %*% v1),
                   as.vector(t(A_current) %*% v2),
                   as.vector(t(A_current) %*% v3),
                   as.vector(t(A_current) %*% v4)
    )
  }


  # All the operator forms involving matrices Q are equivalent to applying the same method for BaseCGF with Q_inner = A^T Q A
  K4operatorAABBfun <- function(tvec, parameter_vector, Q1, Q2) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    cgf$K4operatorAABB(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner)
  }

  K3K3operatorAABBCCfun <- function(tvec, parameter_vector, Q1, Q2, Q3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    # Q3_inner <- tA %*% Q3 %*% A_current
    cgf$K3K3operatorAABBCC(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner, Q1_inner)
  }

  K3K3operatorABCABCfun <- function(tvec, parameter_vector, Q1, Q2, Q3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    Q1_inner <- tA %*% Q1 %*% A_current
    # Q2_inner <- tA %*% Q2 %*% A_current
    # Q3_inner <- tA %*% Q3 %*% A_current
    cgf$K3K3operatorABCABC(as.vector(tA %*% tvec), parameter_vector, Q1_inner, Q1_inner, Q1_inner)
  }

  #### We avoid the factored forms for now (avoiding the potentially expensive loops)
  func_Tfun <- function(tvec, parameter_vector) {
    Q <- solve(K2fun(tvec, parameter_vector))
    K3K3operatorABCABC_val <- K3K3operatorABCABCfun(tvec, parameter_vector, Q, Q, Q)
    K3K3operatorAABBCC_val <- K3K3operatorAABBCCfun(tvec, parameter_vector, Q, Q, Q)
    K4operatorAABB_val <- K4operatorAABBfun(tvec, parameter_vector, Q, Q)
    K4operatorAABB_val/8 - K3K3operatorAABBCC_val/8 - K3K3operatorABCABC_val/12
  }


  # For the factored forms where Q = B D B^T and D has diagonal vector d, note that Q_inner = A^T Q A = (A^T B) D (A^T B)^T
  # Note about sizes: if A is n-by-m then B is n-by-r for some r, and A^T B is m-by-r
  K4operatorAABB_factored <- cgf$.get_private_method("K4operatorAABB_factored")
  K4operatorAABB_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    K4operatorAABB_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2)
  }

  K3K3operatorAABBCC_factored <- cgf$.get_private_method("K3K3operatorAABBCC_factored")
  K3K3operatorAABBCC_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2, B3, d3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    B3_inner <- tA %*% B3
    K3K3operatorAABBCC_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2, B3_inner, d3)
  }

  K3K3operatorABCABC_factored <- cgf$.get_private_method("K3K3operatorABCABC_factored")
  K3K3operatorABCABC_factoredfun <- function(tvec, parameter_vector, B1, d1, B2, d2, B3, d3) {
    A_current <- get_sparse_A(parameter_vector)
    tA <- t(A_current)
    B1_inner <- tA %*% B1
    B2_inner <- tA %*% B2
    B3_inner <- tA %*% B3
    K3K3operatorABCABC_factored(as.vector(tA %*% tvec), parameter_vector, B1_inner, d1, B2_inner, d2, B3_inner, d3)
  }

  # inequality constraints for the transformed variable Y = A * X are the same as those
  # for the original variable X, evaluated at the transformed input A.transpose() * tvec.
  ineq_constraintfun <- function(tvec, parameter_vector) {
    A_current <- get_sparse_A(parameter_vector)
    cgf$ineq_constraint(as.vector(t(A_current) %*% tvec), parameter_vector)
  }


  simulate_fun <- NULL
  if (isTRUE(cgf$has_simulate)) {
    simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      A_current <- get_sparse_A(parameter_vector)
      d_out <- nrow(A_current)
      d_in  <- ncol(A_current)

      # # single-block only
      # if (vector_length != d_out) stop("linearlyMappedCGF$rsim: 'vector_length' must equal nrow(A).", call. = FALSE)
      #

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
  # if (isTRUE(cgf$has_simulate)) {
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

  # ------------------------------------------------------------------
  # # Build the new mapped CGF using createCGF
  # ------------------------------------------------------------------
  createCGF(
    K = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3operatorfun,
    K4operator = K4operatorfun,
    ineq_constraint = ineq_constraintfun,
    analytic_tvec_hat = NULL,
    tilting_exponent = tiltingfun,
    # neg_ll = negllfun,
    func_T = func_Tfun,
    rsim = simulate_fun,
    K4operatorAABB = K4operatorAABBfun,
    K3K3operatorAABBCC = K3K3operatorAABBCCfun,
    K3K3operatorABCABC = K3K3operatorABCABCfun,
    K4operatorAABB_factored = K4operatorAABB_factoredfun,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factoredfun,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factoredfun,
    K2operator = K2operatorfun,
    K2operatorAK2AT = K2operatorAK2ATfun,
    op_name = c(cgf$call_history, "linearlyMappedCGF"),
    ...
  )
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
