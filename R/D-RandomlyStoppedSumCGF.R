# ------------------------------------------------------------------
#  R/D-RandomlyStoppedSumCGF.R
#  Main function: randomlyStoppedSumCGF (exported)
#
#  Purpose
#  -------
#  Create a CGF object for a randomly-stopped sum
#
#      Y = sum_{i=1}^N X_i,
#
#  where:
#    - N is a scalar count random variable with CGF object 'count_cgf'
#    - X_i are i.i.d. copies of a (possibly vector-valued) summand with CGF "summand_cgf"
#    - N is independent of all X_i
#
#  Replication semantics
#  --------------------------------
#  The arguments 'iidReps' and "block_size' refer to replication of the variable Y
#  (i.e. independent observations Y_1,...,Y_B), not replication of summands inside the random sum.
#
#  Because many summand CGFs in this package are vectorized, failing to specify either:
#     - block_size (dimension of one observation Y), or
#     - iidReps (number of observations B),
#  can silently change the statistical model. Therefore this constructor enforces:
#
#     (A) At least one of `block_size` or `iidReps` must be supplied.
#
#     (B) If `iidReps` is NULL (the default), then `block_size` must be supplied and we set iidReps="any".
#         This means: for any tvec whose length is a multiple of block_size, we interpret it as B blocks
#         with B = length(tvec) / block_size.
#
#     (C) If `iidReps` is "any", then `block_size` must be supplied (otherwise we cannot split tvec).
#
#     (D) If `iidReps` is a positive integer, then `block_size` may be NULL (it will be inferred at
#         evaluation as length(tvec) / iidReps), but providing `block_size` explicitly is safer.
#
#  This file also provides efficient implementations of:
#    - logdetK2(t,theta)
#    - K2_solve(t,theta,rhs)
#    - K4operatorAABB(t,theta,Q)
#    - K3K3operatorAABBCC(t,theta,Q)
#    - K3K3operatorABCABC(t,theta,Q)
#    - func_T(t,theta) (uses the fast operator paths)
# ------------------------------------------------------------------








# ------------------------------------------------------------------
# Internal constructor: builds the CGF for ONE replicate of Y.
# Optional iid replication is handled in the exported wrapper.
# ------------------------------------------------------------------
.randomlyStoppedSumCGF_internal <- function(count_cgf, summand_cgf, ...) {


  .trace_mat <- function(M) {
    # Trace = sum of diagonal entries
    sum(diag(M))
  }

  .factor_Q_pd <- function(Q) {
    # For symmetric PD Q, build a factorization
    #   Q = A %*% diag(d) %*% t(A)
    # matching the convention used elsewhere in the package.
    #
    # This is only used for PD inputs.
    U <- chol(Q)                 # upper triangular, Q = t(U) %*% U
    diagU <- diag(U)
                  # # Avoid dividing by 0 if something went very wrong numerically
                  # if (any(diagU == 0)) stop("chol(Q) has a zero diagonal; Q may not be PD.")
    # d <- diagU^2
    d <- diagU * diagU
    A <- t(U) %*% diag(1 / diagU)
    list(A = A, d = d)
  }

  .T3_AAB <- function(cgf, tvec, param, Q, v) {
    # Contraction of the 3rd derivative tensor with a matrix on its first
    # two slots and a vector on the third:
    #
    #   T3_AAB(Q,v) = \sum_{i,j,k} K^{(3)}_{i j k}(t) Q_{i j} v_k.
    #
    # For PD Q, we compute this as a rank-d sum using a Cholesky-based
    # factorization Q = A diag(d) A^T:
    #
    #   T3_AAB(Q,v) = \sum_m d_m K3operator(a_m, a_m, v).
    fac <- .factor_Q_pd(Q)
    A <- fac$A
    d <- fac$d
    out <- 0
    for (m in seq_along(d)) {
      out <- out + d[m] * cgf$K3operator(tvec, param, A[, m], A[, m], v)
    }
    out
  }



  # ------------------------------------------------------------------
  # The CGF composition is:
  #   K_Y(t) = K_N( K_X(t) )
  # where K_X(t) is scalar (even if t is a vector).
  # ------------------------------------------------------------------

  # Convenience: count derivatives at the composed scalar argument
  .count_derivs <- function(s_scalar, param) {
    a1 <- as.numeric(count_cgf$K1(s_scalar, param))[1]
    a2 <- as.numeric(count_cgf$K2(s_scalar, param))[1]
    a3 <- as.numeric(count_cgf$K3operator(s_scalar, param, 1, 1, 1))
    a4 <- as.numeric(count_cgf$K4operator(s_scalar, param, 1, 1, 1, 1))
    list(a1 = a1, a2 = a2, a3 = a3, a4 = a4)
  }

  # -------------------- core vectorized derivatives -----------------

  # K_Y(t) = K_N( K_X(t) )
  K <- function(tvec, param) {
    s <- summand_cgf$K(tvec, param)
    count_cgf$K(s, param)
  }

  # Del K_Y(t) = K_N'(K_X(t)) * ∇K_X(t)
  K1 <- function(tvec, param) {
    s <- summand_cgf$K(tvec, param)
    a1 <- as.numeric(count_cgf$K1(s, param))[1]
    a1 * summand_cgf$K1(tvec, param)
  }

  # Del_2 K_Y(t) = K_N''(K_X(t)) ∇K_X ∇K_X^T + K_N'(K_X(t)) ∇²K_X
  K2 <- function(tvec, param) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    g1 <- summand_cgf$K1(tvec, param)
    g2 <- summand_cgf$K2(tvec, param)
    # Ensure g1 behaves like a column vector for tcrossprod
    g1v <- as.numeric(g1)
    a2 * tcrossprod(g1v) + a1 * g2
  }

  # Bilinear form x^T Del_2 K_Y(t) y
  K2operator <- function(tvec, param, x, y) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    g1 <- as.numeric(summand_cgf$K1(tvec, param))
    # a2*(x^T g1)(y^T g1) + a1*(x^T g2 y)
    a2 * sum(g1 * x) * sum(g1 * y) + a1 * summand_cgf$K2operator(tvec, param, x, y)
  }

  # 3rd-derivative contraction K^{(3)}(w1,w2,w3)
  K3operator <- function(tvec, param, w1, w2, w3) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3

    g1 <- as.numeric(summand_cgf$K1(tvec, param))
    # a1 * K3_X(w1,w2,w3)
    out <- a1 * summand_cgf$K3operator(tvec, param, w1, w2, w3)
    # a2 * [ (w1^T g2 w2)(w3^T g1) + perms ]
    out <- out + a2 * (
      summand_cgf$K2operator(tvec, param, w1, w2) * sum(g1 * w3) +
        summand_cgf$K2operator(tvec, param, w1, w3) * sum(g1 * w2) +
        summand_cgf$K2operator(tvec, param, w2, w3) * sum(g1 * w1)
    )
    # a3 * (w1^T g1)(w2^T g1)(w3^T g1)
    out <- out + a3 * sum(g1 * w1) * sum(g1 * w2) * sum(g1 * w3)
    out
  }

  # 4th-derivative contraction K^{(4)}(w1,w2,w3,w4)
  K4operator <- function(tvec, param, w1, w2, w3, w4) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3
    a4 <- dv$a4

    g1 <- as.numeric(summand_cgf$K1(tvec, param))

    # Partition formula (Faà di Bruno for scalar outer function):
    #   h4 = a1*s4
    #      + a2[ sum_{(3,1)} s3*s1 + sum_{(2,2)} s2*s2 ]
    #      + a3[ sum_{(2,1,1)} s2*s1*s1 ]
    #      + a4[ s1*s1*s1*s1 ]

    s1 <- function(w) sum(g1 * w)
    s2 <- function(u, v) summand_cgf$K2operator(tvec, param, u, v)
    s3 <- function(u, v, w) summand_cgf$K3operator(tvec, param, u, v, w)
    s4 <- function(u, v, w, z) summand_cgf$K4operator(tvec, param, u, v, w, z)

    out <- a1 * s4(w1, w2, w3, w4)

    out <- out + a2 * (
      # (3,1) partitions
      s3(w1, w2, w3) * s1(w4) +
        s3(w1, w2, w4) * s1(w3) +
        s3(w1, w3, w4) * s1(w2) +
        s3(w2, w3, w4) * s1(w1) +
        # (2,2) partitions
        s2(w1, w2) * s2(w3, w4) +
        s2(w1, w3) * s2(w2, w4) +
        s2(w1, w4) * s2(w2, w3)
    )

    out <- out + a3 * (
      # (2,1,1) partitions (6 of them)
      s2(w1, w2) * s1(w3) * s1(w4) +
        s2(w1, w3) * s1(w2) * s1(w4) +
        s2(w1, w4) * s1(w2) * s1(w3) +
        s2(w2, w3) * s1(w1) * s1(w4) +
        s2(w2, w4) * s1(w1) * s1(w3) +
        s2(w3, w4) * s1(w1) * s1(w2)
    )

    out <- out + a4 * s1(w1) * s1(w2) * s1(w3) * s1(w4)
    out
  }

  # -------------------- inequality constraints ----------------------

  ineq_constraint <- function(tvec, param) {
    summand_ineq <- summand_cgf$ineq_constraint(tvec, param)
    s <- summand_cgf$K(tvec, param)
    count_ineq <- count_cgf$ineq_constraint(s, param)

    out_ <- numeric(length(count_ineq) + length(summand_ineq)) * param[1]
    if (length(count_ineq) > 0) out_[1:length(count_ineq)] <- count_ineq
    if (length(summand_ineq) > 0)
      out_[(length(count_ineq) + 1):length(out_)] <- summand_ineq
    out_
  }

  # -------------------- efficient K2 solve + logdet -----------------

  # For Y, the Hessian is
  #   K2_Y = a2 * mu mu^T + a1 * Sigma,
  # a rank-1 update of a scaled Sigma.
  # We exploit Sherman-Morrison + matrix determinant lemma.

  K2_solve <- function(tvec, param, rhs) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    # Use summand's own K2_solve (it may be specialized / sparse)
    Sig_inv_rhs <- summand_cgf$K2_solve(tvec, param, rhs)
    Sig_inv_mu  <- summand_cgf$K2_solve(tvec, param, mu)

    # A = a1*Sigma, so A^{-1} = (1/a1)*Sigma^{-1}
    x  <- Sig_inv_rhs / a1
    u  <- Sig_inv_mu  / a1

    q <- sum(mu * Sig_inv_mu)  # mu^T Sigma^{-1} mu
    denom <- 1 + (a2 / a1) * q

    if (is.matrix(x)) {
      cvec <- as.numeric(crossprod(mu, x)) # length = ncol(x)
      x - u %*% matrix((a2 * cvec) / denom, nrow = 1)
    } else {
      cscal <- sum(mu * x)
      x - u * ((a2 * cscal) / denom)
    }
  }

  logdetK2 <- function(tvec, param) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig_inv_mu <- summand_cgf$K2_solve(tvec, param, mu)
    q <- sum(mu * Sig_inv_mu)  # mu^T Sigma^{-1} mu
    d <- length(mu)
    logdet_Sigma <- as.numeric(summand_cgf$logdetK2(tvec, param))
    as.numeric(d * log(a1) + logdet_Sigma + log1p((a2 / a1) * q))
  }


  # -------------------- efficient operators for func_T --------------

  # K4operatorAABB(t,Q) = \sum_{i,j,k,l} K4_{i j k l} Q_{i j} Q_{k l}.
  # The derived closed form avoids the base-class rank-factor triple loops.

  K4operatorAABB <- function(tvec, param, Q) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3
    a4 <- dv$a4

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig <- summand_cgf$K2(tvec, param)

    v <- Q %*% mu
    q <- as.numeric(crossprod(mu, v))

    tr <- .trace_mat(Sig %*% Q)
    tr12 <- .trace_mat(Sig %*% Q %*% Sig %*% Q)

    # T3 contractions: \sum_{i,j,k} K3_{i j k} Q_{i j} v_k
    T3_Q_v <- .T3_AAB(summand_cgf, tvec, param, Q, v)
    
    # Summand's own K4 AABB contraction
    K4_X <- summand_cgf$K4operatorAABB(tvec, param, Q)

    # Cross quadratic term
    quad12 <- as.numeric(crossprod(v, Sig %*% v))

    a1 * K4_X +
      a2 * (4 * T3_Q_v + tr * tr + 2 * tr12) +
      a3 * (2 * tr * q + 4 * quad12) +
      a4 * (q * q)
  }


  K3K3operatorAABBCC <- function(tvec, param, Q) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig <- summand_cgf$K2(tvec, param)

    v <- Q %*% mu
    qmm <- as.numeric(crossprod(mu, v))

    trSQ <- .trace_mat(Sig %*% Q)
    sv <- Sig %*% v
    uU <- mu * trSQ + 2 * sv
    QuU <- Q %*% uU

    # Required summand terms
    T3_Q_v   <- .T3_AAB(summand_cgf, tvec, param, Q, v)
    T3_Q_QuU <- .T3_AAB(summand_cgf, tvec, param, Q, QuU)

    base_T3T3 <- summand_cgf$K3K3operatorAABBCC(tvec, param, Q)

    UU <- as.numeric(crossprod(uU, Q %*% uU))
    UM <- qmm * as.numeric(crossprod(uU, v))
    MM <- (qmm * qmm) * qmm

    (a1 * a1) * base_T3T3 +
      2 * a1 * a2 * T3_Q_QuU +
      2 * a1 * a3 * qmm * T3_Q_v +
      (a2 * a2) * UU +
      2 * a2 * a3 * UM +
      (a3 * a3) * MM
  }

  K3K3operatorABCABC <- function(tvec, param, Q) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig <- summand_cgf$K2(tvec, param)

    v <- Q %*% mu
    qmm <- as.numeric(crossprod(mu, v))

    # Summand tensor contractions
    base_T3T3 <- summand_cgf$K3K3operatorABCABC(tvec, param, Q)
    T3_vvv    <- summand_cgf$K3operator(tvec, param, v, v, v)

    # Cross term between T3_X and (Sigma,mu) part:
    # needs T3_AAB(M,v) with M = Q Sigma Q (symmetric PD if Q,Sigma PD)
    M <- Q %*% Sig %*% Q
    T3_M_v <- .T3_AAB(summand_cgf, tvec, param, M, v)

    # (Sigma,mu) self term (no Cholesky needed):
    # ||sym(Sigma,mu)||^2_{Q,Q,Q} = 3 tr(Sigma Q Sigma Q) * (mu^T Q mu)
    #                             + 6 (Sigma v)^T Q (Sigma v).
    tr_SQSQ <- .trace_mat(Sig %*% Q %*% Sig %*% Q)
    sv <- Sig %*% v
    sv_Q_sv <- as.numeric(crossprod(sv, Q %*% sv))
    UU <- 3 * tr_SQSQ * qmm + 6 * sv_Q_sv

    # Cross term between (Sigma,mu) and mu^\otimes3:
    quad_vSv <- as.numeric(crossprod(v, Sig %*% v))
    UM <- 3 * quad_vSv * qmm

    MM <- (qmm * qmm) * qmm

    (a1 * a1) * base_T3T3 +
      6 * a1 * a2 * T3_M_v +
      2 * a1 * a3 * T3_vvv +
      (a2 * a2) * UU +
      2 * a2 * a3 * UM +
      (a3 * a3) * MM
  }


  func_T <- function(tvec, param) {
    d <- length(tvec)

    Q <- K2_solve(tvec, param, diag(d))  # Q = K2^{-1}

    K4_AABB <- K4operatorAABB(tvec, param, Q)
    K3K3_AABBCC <- K3K3operatorAABBCC(tvec, param, Q)
    K3K3_ABCABC <- K3K3operatorABCABC(tvec, param, Q)
    (K4_AABB / 8) - (K3K3_AABBCC / 8) - (K3K3_ABCABC / 12)
  } ### Needed? May overlap with default

  combined_history <- paste0(
    "count: ", count_cgf$call_history, "\n",
    "summand: ", summand_cgf$call_history
  )
  op_name_vec <- c(combined_history, "randomlyStoppedSumCGF")











  rsim <- NULL
  if (isTRUE(count_cgf$has_rsim) && isTRUE(summand_cgf$has_rsim)) {
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL,
                     max_total_summands = NULL,
                     ...) {

      d <- as.integer(vector_length)

      # Under tilt t: X_i are tilted by t; N is tilted by s = K_X(t) (scalar)
      count_tvec <- NULL
      if (!is.null(tvec)) {
        s <- summand_cgf$K(tvec, parameter_vector)
        s <- as.numeric(s)
        if (length(s) != 1L || !is.finite(s)) {
          stop("randomlyStoppedSumCGF$rsim: K_X(t) must be a finite scalar.", call. = FALSE)
        }
        count_tvec <- s
      }

      # Sample counts N_1,...,N_n
      N_vec <- count_cgf$rsim(
        n = n,
        vector_length = 1L,
        parameter_vector = parameter_vector,
        tvec = count_tvec,
        flatten = TRUE,
        ...
      )
      N_vec <- as.numeric(N_vec)
      if (any(!is.finite(N_vec))) {
        stop("randomlyStoppedSumCGF$rsim: count_cgf$rsim() returned non-finite values.", call. = FALSE)
      }

      # Enforce integer-ish and >= 0
      tol <- 1e-8
      N_round <- round(N_vec)
      if (any(abs(N_vec - N_round) > tol)) {
        stop(
          "randomlyStoppedSumCGF$rsim: count_cgf$rsim() returned non-integer counts. ",
          "First few draws: ", paste(utils::head(N_vec), collapse = ", "),
          call. = FALSE
        )
      }
      N_int <- as.integer(N_round)
      if (any(N_int < 0L)) {
        stop("randomlyStoppedSumCGF$rsim: counts must be >= 0.", call. = FALSE)
      }

      # total number of summands to generate
      N_total_num <- sum(as.double(N_int))
      if (!is.finite(N_total_num) || N_total_num < 0) {
        stop("randomlyStoppedSumCGF$rsim: invalid total count.", call. = FALSE)
      }

      if (N_total_num > .Machine$integer.max) {
        stop(
          "randomlyStoppedSumCGF$rsim: total number of summands is too large (",
          N_total_num, ").",
          call. = FALSE
        )
      }
      N_total <- as.integer(N_total_num)

      # safety guard against huge allocations
      if (!is.null(max_total_summands)) {
        if (length(max_total_summands) != 1L || !is.finite(max_total_summands) ||
            max_total_summands < 0L || max_total_summands != as.integer(max_total_summands)) {
          stop("'max_total_summands' must be NULL or a nonnegative integer.", call. = FALSE)
        }
        mts <- as.integer(max_total_summands)
        if (N_total > mts) {
          stop(
            "randomlyStoppedSumCGF$rsim: total number of summands (", N_total,
            ") exceeds max_total_summands (", mts, ").",
            call. = FALSE
          )
        }
      }

      # Edge-case: all sampled counts are zero => Y is identically 0
      if (N_total == 0L) return(matrix(0, nrow = d, ncol = n))

      # sample all summands at once
      X_all <- summand_cgf$rsim(
        n = N_total,
        vector_length = d,
        parameter_vector = parameter_vector,
        tvec = tvec,
        flatten = FALSE,
        ...
      )

      # Segment-sum into n outputs
      Y <- matrix(0, nrow = d, ncol = n)
      pos <- 1L
      for (j in seq_len(n)) {
        nj <- N_int[j]
        if (nj > 0L) {
          if (nj == 1L) {
            Y[, j] <- X_all[, pos]
          } else {
            Y[, j] <- rowSums(X_all[, pos:(pos + nj - 1L), drop = FALSE])
          }
          pos <- pos + nj
        }
      }
      Y
    }
  }

  # Build args list (names match createCGF parameters exactly)
  cgf_args <- list(
    K = K,
    K1 = K1,
    K2 = K2,
    K2operator = K2operator,
    K3operator = K3operator,
    K4operator = K4operator,
    K2_solve = K2_solve,
    logdetK2 = logdetK2,
    rsim = rsim,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    func_T = func_T,
    ineq_constraint = ineq_constraint,
    op_name = op_name_vec
  )

  do.call(createCGF, c(cgf_args, list(...)))
}










#' @title CGF for a randomly-stopped sum
#'
#' @description
#' Builds a CGF object for the random vector
#' \deqn{Y = \sum_{i=1}^{N} X_i,}
#' where:
#' \itemize{
#'   \item \eqn{N} is a non-negative integer-valued scalar random variable with CGF \code{count_cgf},
#'   \item \eqn{X_i} are i.i.d. copies of a (possibly vector-valued) summand with CGF \code{summand_cgf},
#'   \item \eqn{N} is independent of all \eqn{X_i}.
#' }
#'
#'
#' @param count_cgf A \code{CGF} object for the (scalar) count variable \eqn{N}.
#' @param summand_cgf A \code{CGF} object for the summand \eqn{X}.
#' @param block_size Optional. The dimension \eqn{d} of one observation of \eqn{Y} (and hence of \eqn{X}).
#'   This can be a positive integer (fixed \eqn{d}) or \code{NULL}
#' @param iidReps Optional. Replication count \eqn{B} for i.i.d. observations \eqn{Y_1,\ldots,Y_B}.
#'   May be a positive integer, \code{"any"}, or \code{NULL}.
#' @param ... Additional named arguments passed to CGF creation.
#'
#'
#'
#' @details
#' **What do `block_size` and `iidReps` mean for RSS models?**
#'
#' The replication arguments refer to i.i.d. replication of the variable \eqn{Y},
#' not replication of summands inside the random sum.
#'
#' If you have \eqn{B} independent observations,
#' \deqn{Y_b = \sum_{i=1}^{N_b} X_{b,i}, \qquad b=1,\ldots,B,}
#' then:
#' \itemize{
#'   \item \code{block_size} is the dimension of one \eqn{Y_b} (which equals the dimension of one \eqn{X_{b,i}}).
#'   \item \code{iidReps} is the number of independent observations \eqn{B}.
#' }
#'
#' **Why do we require `block_size` and/or `iidReps`?**
#'
#' Many summand CGFs in this package are vectorized: if you accidentally feed a long vector \code{tvec}
#' into an 'untagged' RSS CGF, the summand CGF may interpret that as a higher-dimensional summand \eqn{X}
#' (or multiple i.i.d. summands inside one \eqn{X}), which silently changes the model.
#'
#' To prevent this, \code{randomlyStoppedSumCGF()} requires at least one of \code{block_size} or \code{iidReps}.
#'
#' The semantics are:
#' \itemize{
#'   \item If \code{iidReps} is \code{NULL} (default), then \code{block_size} must be provided and we set \code{iidReps="any"}.
#'         This means the number of observations is inferred at evaluation time as \code{B = length(tvec)/block_size}.
#'   \item If \code{iidReps} is \code{"any"}, then \code{block_size} must be provided.
#'   \item If \code{iidReps} is a positive integer, then \code{block_size} may be \code{NULL} (it will be inferred as
#'         \code{length(tvec)/iidReps}), but providing \code{block_size} explicitly is encouraged.
#' }
#'
#' @examples
#' \dontrun{
#' ## Scalar RSS, B observations:
#' K.N <- GeometricModelCGF(prob = adaptor(fixed_param = 0.4))
#' K.X <- BinomialModelCGF(n = adaptor(fixed_param = 1), p = adaptor(indices = 1))
#'
#' # Either specify B explicitly...
#' K.U1 <- randomlyStoppedSumCGF(K.N, K.X, iidReps = 40)      # block_size inferred as length(tvec)/40
#'
#' # ...or specify the scalar block size and let B be inferred from length(tvec):
#' K.U2 <- randomlyStoppedSumCGF(K.N, K.X, block_size = 1)    # iidReps defaults to "any"
#'
#' ## Vector RSS: summand is d-dimensional
#' d <- 3
#' K.Xv <- MultinomialModelCGF(n = adaptor(fixed_param = 1),
#'                             prob_vec = adaptor(indices = 1:d))
#' K.Y  <- randomlyStoppedSumCGF(PoissonCGF, K.Xv, block_size = d, iidReps = "any")
#' }
#'
#' @return A `CGF` object.
#' @export
randomlyStoppedSumCGF <- function(count_cgf,
                                  summand_cgf,
                                  block_size = NULL,
                                  iidReps = NULL,
                                  ...) {
  if (!inherits(count_cgf, "CGF")) stop("'count_cgf' must be a CGF object.")
  if (!inherits(summand_cgf, "CGF")) stop("'summand_cgf' must be a CGF object.")
  if (!is.null(iidReps)) .check_iidReps(iidReps)

  # out_cgf <- .randomlyStoppedSumCGF_internal(count_cgf, summand_cgf, ...)
  out_cgf <- .randomlyStoppedSumCGF_internal(
    count_cgf, summand_cgf,
    ...
  )


  # require at least one replication hint to avoid silent mistakes.
  if (is.null(block_size) && is.null(iidReps)) {
    stop("randomlyStoppedSumCGF(): Please supply at least one of 'block_size' or 'iidReps' to disambiguate i.i.d. replication semantics.")
  }

  if (!is.null(iidReps) && iidReps == 1 && is.null(block_size)) return(out_cgf)

  # if (!is.null(iidReps) && iidReps == 1 && is.null(block_size)) {
  #   # Allow explicit iidReps=1 as a way to declare "single replicate".
  #   return(out_cgf)
  # }

  if (is.null(iidReps)) iidReps <- "any"

  # If iidReps="any", block_size is required
  if (identical(iidReps, "any") && is.null(block_size)) {
    stop("randomlyStoppedSumCGF(): iidReps='any' requires a non-NULL 'block_size' ")
  }

  iidReplicatesCGF(cgf = out_cgf, iidReps = iidReps, block_size = block_size)
}
