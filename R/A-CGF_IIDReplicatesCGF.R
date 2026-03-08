# R/IIDReplicatesCGF.R
# Objects: iidReplicatesCGF



.iidReplicatesCGF_internal <- function(cgf, iidReps, block_size) {

  # d = block_size, B = number of blocks, N = length of tvec
  chunkIndices <- function(i, block_size) {
    seq.int((i - 1)*block_size + 1, i*block_size)
  }

  # format block_size for labels; never forces evaluation
  .format_block_size_tag <- function(block_size) {
    if (is.null(block_size)) return(NULL)
    if (is.function(block_size)) {
      lab <- attr(block_size, "label")
      if (is.character(lab) && nzchar(lab)) return(lab)
      return("d(theta)")   # default symbolic tag
    }
    block_size
  }



  # fetch some private methods from the base CGF
  base_tilting_exponent <- cgf$.private_api$tilting_exponent
  base_neg_ll <- cgf$.private_api$neg_ll
  base_func_T <- cgf$.private_api$func_T

  # ------------------------------------------------------------------
  # Now all methods in a unified manner
  # ------------------------------------------------------------------


  # K => sum over blocks
  K <- function(tvec, param) {
    N <- length(tvec)

    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])

    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      # s <- idxs[i,1]; e <- idxs[i,2]
      total <- total + cgf$K(tvec[idx], param)
    }
    total
  }


  # ------------------------------------------------------------------
  # #### NOTE:
  # The line `out_ <- numeric(n)` currently triggers an error during
  # Tape recording. To address this issue, there are two potential solutions:
  #
  # 1. Utilize the ADoverload package to overload the `[<-` operator.
  # 2. Modify the line to create an explicit dependency on a parameter,
  #    such as `out_ <- numeric(n) * param[1]`.
  #
  # For the time being, we use the second approach to avoid any
  # complications with operator overloading.
  # ------------------------------------------------------------------
  # K1 => piecewise concatenation
  K1 <- function(tvec, param) {
    N <- length(tvec)

    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    out_ <- numeric(N) * param[1]
    for (i in seq_len(B)) {
      # s <- idxs[i,1]; e <- idxs[i,2]
      # out[s:e] <- cgf$K1(tvec[s:e], param)
      idx <- chunkIndices(i, d)
      out_[idx] <- cgf$K1(tvec[idx], param)
    }
    out_
  }


  # K1fun <- function(tvec, param) {
  #   n <- length(tvec)
  #   N <- get_nBlocks(n)
  #   bS <- get_blockSize(n)
  #   out_ <- numeric(n) * param[1]  # tie to 'param' to avoid tape error
  #   for (i in seq_len(N)) {
  #     idx <- chunkIndices(i, bS)
  #     out_[idx] <- cgf$K1(tvec[idx], param)
  #   }
  #   out_
  # }



  # K2 => block-diagonal
  K2 <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    accum <- matrix(0, nrow = N, ncol = N) * param[1]
    for (i in seq_len(B)) {
      # s <- idxs[i, 1]; e <- idxs[i, 2]
      idx <- chunkIndices(i, d)
      k2 <- cgf$K2(tvec[idx], param)
      if (is.null(dim(k2))) k2 <- matrix(k2, nrow = d, ncol = d)  # d = 1 hardening
      accum[idx, idx] <- as.matrix(k2)
      # s <- idxs[i,1]; e <- idxs[i,2]
      # accum[s:e, s:e] <- cgf$K2(tvec[s:e], param)
    }
    accum
  }




  # tilting_exponent => sum
  tilting_exponent <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + base_tilting_exponent(tvec[idx], param)
    }
    total
  }


  # neg_ll => sum
  neg_ll <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + base_neg_ll(tvec[idx], param)
    }
    total
  }


  # func_T => sum
  func_T <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + base_func_T(tvec[idx], param)
    }
    total
  }

  # K2operator => sum
  K2operator <- function(tvec, param, x, y) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + cgf$K2operator(tvec[idx], param, x[idx], y[idx])
    }
    total
  }



  K2operatorAK2AT <- function(tvec, param, Bmat) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    lay <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])

    if (ncol(Bmat) != N) {
      stop("K2operatorAK2AT: Bmat must have ncol == length(tvec). ",
           "Got ncol(Bmat)=", ncol(Bmat), ", length(tvec)=", N, ".")
    }

    r <- nrow(Bmat)
    out <- matrix(0, nrow = r, ncol = r) * param[1]  # keep AD type if needed

    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      out <- out + cgf$K2operatorAK2AT(
        tvec[idx], param,
        Bmat[, idx, drop = FALSE]
      )
    }
    out
  }


  # K3operator => sum
  K3operator <- function(tvec, param, v1, v2, v3) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + cgf$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
    }
    total
  }

  # K4operator => sum
  K4operator <- function(tvec, param, v1, v2, v3, v4) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + cgf$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
    }
    total
  }


  # K4operatorAABB => sum
  K4operatorAABB <- function(tvec, param, Q) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      # slice out sub-block of Q
      ### (Verify?) Q is block-diagonal of size (n·iidReps) × (n·iidReps). The relevant sub-block is n×n
      Qsub <- Q[idx, idx, drop = FALSE]
      total <- total + cgf$K4operatorAABB(tvec[idx], param, Qsub)
    }
    total
  }


  # K3K3operatorAABBCC => sum
  K3K3operatorAABBCC <- function(tvec, param, Q) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      Qsub <- Q[idx, idx, drop = FALSE]
      total <- total + cgf$K3K3operatorAABBCC(tvec[idx], param, Qsub, Qsub, Qsub)
      ############ Note: check, this seems to be incorrect, as Q may have non-zero off-diagonal blocks (cf. comment "Verify?" in K4operatorAABB)
      ##### cf. code for vectorized CGFs
    }
    total
  }

  # K3K3operatorABCABC => sum
  K3K3operatorABCABC <- function(tvec, param, Q) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      Qsub <- 1[idx, idx, drop = FALSE]
      total <- total + cgf$K3K3operatorABCABC(tvec[idx], param, Qsub, Qsub, Qsub)
      ############ Note: check, this seems to be incorrect, as Q may have non-zero off-diagonal blocks (cf. comment "Verify?" in K4operatorAABB)
      ##### cf. code for vectorized CGFs
    }
    total
  }


  # ineq_constraint => concatenation
  ineq_constraint <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
    # d_cur <- .block_size_value(block_size, getval_param)
    lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])



    # call once for length
    idx_first <- seq.int(1, d)
    first_val <- cgf$ineq_constraint(tvec[idx_first], param)
    L <- length(first_val)
    if (L == 0L) return(numeric(0) * param[1])

    out_ <- numeric(L * B) * param[1] #### Modified to depend on param
    out_[1:L] <- first_val
    if (B > 1) {
      for (i in 2:B) {
        idx <- chunkIndices(i, d)
        val <- cgf$ineq_constraint(tvec[idx], param)
        start_ <- (i - 1)*L + 1
        out_[start_:(i*L)] <- val
      }
    }
    out_
  }

  # for the analytic_tvec_hat:
  # We'll do a chunk approach if cgf$analytic_tvec_hat() is non-NULL:
  # e.g. chunk x => pass each chunk to cgf$analytic_tvec_hat => combine?
  #### Check if this doesn't make sense, (default to NULL if that's the case)
  analytic_tvec_hat <- NULL # If the base CGF had no valid function, just return NULL
  if (isTRUE(cgf$has_analytic_tvec_hat)) {
    analytic_tvec_hat <- function(x, param) {
      N <- length(x)
      d_cur <- .block_size_value(block_size, param)
      # getval_param <- ifelse(is(param, "advector"), RTMB:::getValues(param), param)
      # d_cur <- .block_size_value(block_size, getval_param)
      lay  <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
      d    <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])


      out_ <- numeric(N) * param[1]
      for (i in seq_len(B)) {
        idx <- chunkIndices(i, d)
        out_[idx] <- cgf$analytic_tvec_hat(x[idx], param)
      }
      out_
    }
  }

  bs_tag <- .format_block_size_tag(block_size)

  pieces <- character(0)
  if (!identical(iidReps, "any")) pieces <- c(pieces, sprintf("iidReps=%d", as.integer(iidReps)))
  if (!is.null(bs_tag))          pieces <- c(pieces, sprintf("bS=%s", bs_tag))

  op_label <- if (length(pieces)) {
    sprintf("iidReplicatesCGF(%s)", paste(pieces, collapse = ","))
  } else {
    "iidReplicatesCGF"
  }





  K2_solve <- function(tvec, param, rhs) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    lay <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])

    # vector RHS
    if (is.null(dim(rhs))) {
      if (length(rhs) != N) stop("K2_solve: rhs length mismatch.")
      out <- numeric(N) * param[1]   # RTMB tape-safe pattern you already use
      for (i in seq_len(B)) {
        idx <- chunkIndices(i, d)
        out[idx] <- cgf$K2_solve(tvec[idx], param, rhs[idx])
      }
      return(out)
    }

    # matrix RHS
    if (nrow(rhs) != N) stop("K2_solve: rhs nrow mismatch.")
    k <- ncol(rhs)
    out <- matrix(0, nrow = N, ncol = k) * param[1]
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      out[idx, ] <- cgf$K2_solve(tvec[idx], param, rhs[idx, , drop = FALSE])
    }
    out
  }


  logdetK2 <- function(tvec, param) {
    N <- length(tvec)
    d_cur <- .block_size_value(block_size, param)
    lay <- .resolve_rep_layout(N, block_size = d_cur, iidReps = iidReps)
    d <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])

    total <- 0
    for (i in seq_len(B)) {
      idx <- chunkIndices(i, d)
      total <- total + cgf$logdetK2(tvec[idx], param)
    }
    total
  }

  # ------------------------------------------------------------------
  # Optional simulator: forward if the base CGF can simulate.
  # iidReplicatesCGF only changes how long tvec vectors are interpreted (block sums),
  # so simulation can be forwarded directly.
  # ------------------------------------------------------------------
  rsim <- NULL
  if (isTRUE(cgf$has_rsim)) {
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      d_cur <- .block_size_value(block_size, parameter_vector)
      lay <- .resolve_rep_layout(vector_length, block_size = d_cur, iidReps = iidReps)
      d <- as.integer(lay[["d"]]); B <- as.integer(lay[["B"]])

      if (!is.null(tvec) && length(tvec) != vector_length) stop("iidReplicatesCGF$rsim: if supplied, 'tvec' must have length == vector_length.", call. = FALSE)

      out <- matrix(0, nrow = vector_length, ncol = n)

      if (is.null(tvec)) {
        X <- cgf$rsim(
          n = n * B,
          vector_length = d,
          parameter_vector = parameter_vector,
          tvec = NULL,
          flatten = FALSE,
          ...
        )
        for (j in seq_len(n)) {
          cols <- ((j - 1L) * B + 1L):(j * B)
          out[, j] <- as.vector(X[, cols, drop = FALSE])
        }
        return(out)
      }

      for (b in seq_len(B)) {
        # idx <- chunkIndices(b, d)
        idx <- ((b - 1L) * d + 1L):(b * d)
        Xb <- cgf$rsim(
          n = n,
          vector_length = d,
          parameter_vector = parameter_vector,
          tvec = tvec[idx],
          flatten = FALSE,
          ...
        )
        out[idx, ] <- Xb
      }

      out
    }
  }




  # ------------------------------------------------------------------
  # Build the new CGF object
  # ------------------------------------------------------------------
  op_name <- c(cgf$call_history, op_label)

  cgf_args <- list(
    K = K,
    K1 = K1,
    K2 = K2,
    K3operator = K3operator,
    K4operator = K4operator,
    tilting_exponent = tilting_exponent,
    neg_ll = neg_ll,
    func_T = func_T,
    ineq_constraint = ineq_constraint,
    analytic_tvec_hat = analytic_tvec_hat,
    K2operator = K2operator,
    K2operatorAK2AT = K2operatorAK2AT,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K2_solve = K2_solve,
    logdetK2 = logdetK2,
    rsim = rsim,
    op_name = op_name
  )

  do.call(createCGF, cgf_args)

}



















#' @title Replicate a CGF object over multiple i.i.d. blocks and/or with a fixed block size
#'
#' @description
#' Extends a given `CGF` object to handle multiple i.i.d. blocks. You can specify:
#'
#' - \code{iidReps} only,
#' - \code{block_size} only,
#' - or both \code{iidReps} and \code{block_size}.
#'
#' @param cgf A `CGF` object.
#' @param iidReps Either `"any"` (default) or a positive integer.
#' @param block_size Either NULL, a positive integer, or a function
#'   \code{function(param) -> positive integer}. When a function is provided,
#'   the block size is determined dynamically from the current parameter vector.
#'
#' @return A `CGF` object
#'
#' @details
#' Let \eqn{N = \mathrm{length(tvec)}}, \code{d = block_size}, and \code{iidReps} the number of blocks.
#' \itemize{
#'   \item If `iidReps = "any"` and `block_size = NULL`: pass-through (no splitting).
#'   \item If `iidReps = "any"` and `block_size = d`: require `N %% d == 0`, set `iidReps = N/d`.
#'   \item If `iidReps = m` and `block_size = NULL`: require `N %% m == 0`, set `d = N/m`, `iidReps = m`.
#'   \item If `iidReps = m` and `block_size = d`: require `N == d * m`, set `iidReps = m`.
#' }
#'
#' @examples
#' ## Base CGF: univariate Poisson with lambda(theta) = theta[1]
#' pois <- PoissonModelCGF(lambda = function(th) th[1], iidReps = "any")
#' theta <- c(2)  # rate
#'
#' ## iidReps = "any", block_size = NULL
#' pass <- iidReplicatesCGF(pois, iidReps = "any", block_size = NULL)
#' identical(pass, pois)  # TRUE
#' pass$K1(0, theta)  # same as pois$K1(0, theta)
#'
#' ## block_size only
#' tvec <- c(0.1, -0.2, 0.0, 0.3)    # N = 4
#' agg_b <- iidReplicatesCGF(pois, block_size = 2)  # 2 blocks each of length 2
#' k1_b  <- agg_b$K1(tvec, theta)
#' k1_ref <- c(pois$K1(tvec[1:2], theta), pois$K1(tvec[3:4], theta))
#' all.equal(k1_b, k1_ref) # TRUE
#'
#' ## iidReps only: split tvec into exactly iidReps blocks (block size inferred)
#' t6 <- c(0.1, -0.2, 0.0, 0.3, 0.2, -0.1)  # N = 6
#' agg_m <- iidReplicatesCGF(pois, iidReps = 3)     # B = 3, d = N/B = 2
#' K2_m  <- agg_m$K2(t6, theta)
#' # Off-block covariances are zero (block diagonal)
#' is_zero <- function(M) max(abs(M)) < 1e-12
#' is_zero(K2_m[1:2, 3:4]) && is_zero(K2_m[1:2, 5:6]) && is_zero(K2_m[3:4, 5:6])
#'
#' # Both iidReps and block_size: require N == d * B
#' agg_fix <- iidReplicatesCGF(pois, iidReps = 3, block_size = 2)  # N must be 6 here
#' all.equal(agg_fix$K1(t6, theta), agg_m$K1(t6, theta))  # same objects agg_fix/agg_m
#' \dontrun{
#'   # Mismatch example (errors at evaluation, not at construction):
#'   bad <- iidReplicatesCGF(pois, iidReps = 3, block_size = 4)
#'   bad$K1(t6, theta)  # length(t6) = 6 != 3 * 4  ==> error
#' }
#'
#'
#' ## Here the inner Poisson builder requires EXACTLY 2 inner replicates.
#' inner2 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 2)
#' t6 <- c(0.1, -0.2, 0.0, 0.3, 0.2, -0.1)  # N = 6
#' # With block_size = 2, each outer block (length d = 2) satisfies the inner rule.
#' outer <- iidReplicatesCGF(inner2, block_size = 2)  # B = 3 blocks; d = 2 per block
#' outer$K1(t6, theta)  # valid; each block passes inner iidReps = 2 check
#'
#' ## Non-identical setup
#' ## lambda(theta) returns a vector of length 2; inner builder set to "any".
#' nonid <- PoissonModelCGF(lambda = function(th) c(th[1], 3*th[1]), iidReps = "any")
#' t12 <- rep(c(0.05, -0.10, 0.00, 0.20), 3)  # N = 12
#' # With iidReps = 3 (outer), we have 3 blocks each of size 4.
#' # The inner builder sees d = 4 with 2 inner replicates per block.
#' agg_nonid <- iidReplicatesCGF(nonid, iidReps = 3)
#' agg_nonid$K1(t12, theta)
#'
#' @export
iidReplicatesCGF <- function(cgf, iidReps = "any", block_size = NULL) {
  if (!inherits(cgf, "CGF")) stop("'cgf' must be an object of class 'CGF'.")

  .check_iidReps(iidReps); .check_block_size(block_size)


  # If there's no replication to enforce
  if (identical(iidReps, "any") && is.null(block_size)) return(cgf)
  if (is.numeric(iidReps) && iidReps == 1L && is.null(block_size)) return(cgf)


  .iidReplicatesCGF_internal(
    cgf        = cgf,
    iidReps    = iidReps,
    block_size = block_size
  )

}

