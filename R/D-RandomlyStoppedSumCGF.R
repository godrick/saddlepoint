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

  summand_K2_factor <- .K2_factor_method(summand_cgf)
  extra_args <- list(...)
  extra_args <- extra_args[!vapply(extra_args, is.null, logical(1))]
  extra_names <- names(extra_args)
  if (is.null(extra_names)) extra_names <- rep("", length(extra_args))

  summand_K4AABB_factored <-
    summand_cgf$.private_api$K4operatorAABB_factored
  summand_K3K3AABBCC_factored <-
    summand_cgf$.private_api$K3K3operatorAABBCC_factored
  summand_K3K3ABCABC_factored <-
    summand_cgf$.private_api$K3K3operatorABCABC_factored
  summand_K4AABB_safe <-
    .factored_delegate_is_safe(summand_K4AABB_factored)
  summand_K3K3AABBCC_safe <-
    .factored_delegate_is_safe(summand_K3K3AABBCC_factored)
  summand_K3K3ABCABC_safe <-
    .factored_delegate_is_safe(summand_K3K3ABCABC_factored)

  .validate_factored_Q <- function(tvec, B, dvec, where) {
    B_dim <- dim(B)
    if (length(B_dim) != 2L || B_dim[1L] != length(tvec)) {
      stop(where, ": B must have nrow(B) == length(tvec).", call. = FALSE)
    }
    if (B_dim[2L] != length(dvec)) {
      stop(where, ": ncol(B) must equal length(dvec).", call. = FALSE)
    }
    length(dvec)
  }


  .trace_mat <- function(M) {
    # Trace = sum of diagonal entries
    sum(diag(M))
  }

  .T3_AAB <- function(cgf, tvec, param, Q, v) {
    # Contraction of the 3rd derivative tensor with a matrix on its first
    # two slots and a vector on the third:
    #
    #   T3_AAB(Q,v) = \sum_{i,j,k} K^{(3)}_{i j k}(t) Q_{i j} v_k.
    #
    # Expanding one index in the coordinate basis is exact for every Q:
    #
    #   T3_AAB(Q,v) = \sum_j K3operator(Q[,j], e_j, v).
    #
    # This also covers a PSD singular Q, which occurs when a singular child
    # covariance is rank-completed by the count-variance term.
    d <- nrow(Q)
    basis <- diag(d)
    out <- .ad_zero_scalar(param)
    for (j in seq_len(d)) {
      out <- out + cgf$K3operator(tvec, param, Q[, j], basis[, j], v)
    }
    out
  }

  .T3_AAB_factored <- function(cgf, tvec, param, B, dvec, v) {
    out <- .ad_zero_scalar(param)
    for (j in seq_along(dvec)) {
      bj <- as.vector(B[, j])
      out <- out + dvec[j] * cgf$K3operator(tvec, param, bj, bj, v)
    }
    out
  }

  .coordinate_child_AABBCC <- function(tvec, param, Q) {
    d <- length(tvec)
    slices <- .extract_K3_slices(
      summand_cgf$K3operator, tvec, param, d, diag(1, d)
    )
    u <- .k3_slices_to_aabbcc_vector(slices, Q, param)
    sum(u * as.vector(Q %*% u))
  }

  .coordinate_child_K4_AABB <- function(tvec, param, Q) {
    d <- nrow(Q)
    basis <- diag(d)
    out <- .ad_zero_scalar(param)
    for (i in seq_len(d)) {
      for (j in seq_len(d)) {
        out <- out + summand_cgf$K4operator(
          tvec, param, Q[, i], basis[, i], Q[, j], basis[, j]
        )
      }
    }
    out
  }

  .coordinate_child_ABCABC <- function(tvec, param, Q) {
    d <- length(tvec)
    slices <- .extract_K3_slices(
      summand_cgf$K3operator, tvec, param, d, diag(1, d)
    )
    .k3_slices_abcabc_from_dense_Q(
      list(slices), list(seq_len(d)), Q, param
    )
  }

  # Contract T3_X with M = Q Sigma Q and v without factoring M.  For a
  # genuinely thin Q = B D B', expand in factor space.  At large rank the
  # coordinate expansion uses fewer K3 calls and keeps the established dense
  # arithmetic complexity; M may be singular because .T3_AAB never calls chol.
  .T3_QSigmaQ_v_factored <- function(cgf, tvec, param, B, dvec, G, v) {
    r <- length(dvec)
    out <- .ad_zero_scalar(param)
    if (r == 0L) return(out)

    if (as.double(r) * r <= nrow(B)) {
      for (i in seq_len(r)) {
        bi <- as.vector(B[, i])
        out <- out + ((dvec[i] * G[i, i]) * dvec[i]) *
          cgf$K3operator(tvec, param, bi, bi, v)
        if (i < r) {
          for (j in seq.int(i + 1L, r)) {
            bj <- as.vector(B[, j])
            coefficient <-
              (dvec[i] * G[i, j]) * dvec[j] +
              (dvec[i] * G[j, i]) * dvec[j]
            out <- out + coefficient *
              cgf$K3operator(tvec, param, bi, bj, v)
          }
        }
      }
      return(out)
    }

    M <- B %*% (dvec * (G %*% (dvec * t(B))))
    .T3_AAB(cgf, tvec, param, M, v)
  }

  .project_summand_K2 <- function(tvec, param, B) {
    if (is.null(summand_K2_factor)) {
      return(summand_cgf$K2operatorAK2AT(tvec, param, t(B)))
    }

    r <- ncol(B)
    terms <- summand_K2_factor(tvec, param, t(B))
    if (!is.list(terms) || length(terms) == 0L) {
      stop("Summand K2 factorization supplied no covariance terms.", call. = FALSE)
    }

    G <- .ad_zero_array(c(r, r), param)
    for (term in terms) {
      if (!is.null(term$S)) {
        if (!identical(dim(term$S), c(r, r))) {
          stop("Invalid summand K2 factor term dimensions.", call. = FALSE)
        }
        G <- G + term$S
      } else if (!is.null(term$B) && !is.null(term$d)) {
        term_d <- as.vector(term$d)
        if (nrow(term$B) != r || ncol(term$B) != length(term_d)) {
          stop("Invalid summand K2 factor term dimensions.", call. = FALSE)
        }
        G <- G + term$B %*% (term_d * t(term$B))
      } else {
        stop("Invalid summand K2 factor term.", call. = FALSE)
      }
    }
    G
  }

  # Quantities shared by the three RSS contractions when Q = B D B'.
  # G is only r-by-r.  A child's factor capability is projected directly;
  # otherwise its public covariance sandwich remains authoritative.
  .factored_contraction_terms <- function(tvec, param, B, dvec) {
    s <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    mu <- as.vector(summand_cgf$K1(tvec, param))
    B_mu <- as.vector(t(B) %*% mu)
    weighted_B_mu <- dvec * B_mu
    v <- as.vector(B %*% weighted_B_mu)
    G <- .project_summand_K2(tvec, param, B)
    if (inherits(G, "denseMatrix") && !inherits(G, "adsparse")) {
      G <- as.matrix(G)
    }
    G_weighted_B_mu <- as.vector(G %*% weighted_B_mu)
    weighted_G <- dvec * G

    list(
      dv = dv,
      B_mu = B_mu,
      weighted_B_mu = weighted_B_mu,
      v = v,
      G = G,
      G_weighted_B_mu = G_weighted_B_mu,
      qmm = sum(B_mu * weighted_B_mu),
      trSQ = sum(dvec * diag(G)),
      tr_SQSQ = sum(weighted_G * t(weighted_G)),
      quad_vSv = sum(weighted_B_mu * G_weighted_B_mu)
    )
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

  K2_factor <- NULL
  if (!is.null(summand_K2_factor)) {
    K2_factor <- function(tvec, param, B) {
      s <- summand_cgf$K(tvec, param)
      dv <- .count_derivs(s, param)
      terms <- .K2_factor_scale(
        summand_K2_factor(tvec, param, B),
        dv$a1
      )
      mu <- as.vector(B %*% summand_cgf$K1(tvec, param))
      c(terms, .K2_factor_term(matrix(mu, ncol = 1L), dv$a2))
    }
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

  K2_solve <- function(tvec, param, rhs) {
    if (!is.null(K2_factor)) {
      return(.K2_factor_solve(K2_factor, tvec, param, rhs))
    }

    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    # Preserve a child's specialized / sparse solve when no factor capability
    # is available.  This is the original Sherman-Morrison fast route.
    Sig_inv_rhs <- summand_cgf$K2_solve(tvec, param, rhs)
    Sig_inv_mu  <- summand_cgf$K2_solve(tvec, param, mu)

    x <- Sig_inv_rhs / a1
    u <- Sig_inv_mu / a1
    q <- sum(mu * Sig_inv_mu)
    denom <- 1 + (a2 / a1) * q

    if (is.matrix(x)) {
      cvec <- as.numeric(crossprod(mu, x))
      x - u %*% matrix((a2 * cvec) / denom, nrow = 1)
    } else {
      cscal <- sum(mu * x)
      x - u * ((a2 * cscal) / denom)
    }
  }

  logdetK2 <- function(tvec, param) {
    if (!is.null(K2_factor)) {
      return(.K2_factor_logdet(K2_factor, tvec, param))
    }

    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig_inv_mu <- summand_cgf$K2_solve(tvec, param, mu)
    q <- sum(mu * Sig_inv_mu)
    d <- length(mu)
    logdet_Sigma <- as.numeric(summand_cgf$logdetK2(tvec, param))
    as.numeric(d * log(a1) + logdet_Sigma + log1p((a2 / a1) * q))
  }


  # -------------------- efficient operators for func_T --------------

  # K4operatorAABB(t,Q) = \sum_{i,j,k,l} K4_{i j k l} Q_{i j} Q_{k l}.
  # The derived closed form avoids the base-class rank-factor triple loops.

  K4operatorAABB_core <- function(tvec, param, Q, K4_X) {
    s  <- summand_cgf$K(tvec, param)
    dv <- .count_derivs(s, param)
    a1 <- dv$a1
    a2 <- dv$a2
    a3 <- dv$a3
    a4 <- dv$a4

    mu <- as.numeric(summand_cgf$K1(tvec, param))
    Sig <- summand_cgf$K2(tvec, param)

    v <- Q %*% mu
    qmm <- as.numeric(crossprod(mu, v))

    trSQ <- .trace_mat(Sig %*% Q)
    tr_SQSQ <- .trace_mat(Sig %*% Q %*% Sig %*% Q)

    # T3 contractions: \sum_{i,j,k} K3_{i j k} Q_{i j} v_k
    T3_Q_v <- .T3_AAB(summand_cgf, tvec, param, Q, v)
    # Cross quadratic term
    quad_vSv <- as.numeric(crossprod(v, Sig %*% v))

    a1 * K4_X +
      a2 * (4 * T3_Q_v + trSQ * trSQ + 2 * tr_SQSQ) +
      a3 * (2 * trSQ * qmm + 4 * quad_vSv) +
      a4 * (qmm * qmm)
  }

  K4operatorAABB <- function(tvec, param, Q) {
    K4operatorAABB_core(
      tvec,
      param,
      Q,
      summand_cgf$K4operatorAABB(tvec, param, Q)
    )
  }

  K3K3operatorAABBCC_core <- function(tvec, param, Q, base_T3T3) {
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

  K3K3operatorAABBCC <- function(tvec, param, Q) {
    K3K3operatorAABBCC_core(
      tvec,
      param,
      Q,
      summand_cgf$K3K3operatorAABBCC(tvec, param, Q)
    )
  }

  K3K3operatorABCABC_core <- function(tvec, param, Q, base_T3T3) {
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
    T3_vvv    <- summand_cgf$K3operator(tvec, param, v, v, v)

    # Cross term between T3_X and (Sigma,mu) part. M can be singular even
    # when the final RSS covariance is positive definite.
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

  K3K3operatorABCABC <- function(tvec, param, Q) {
    K3K3operatorABCABC_core(
      tvec,
      param,
      Q,
      summand_cgf$K3K3operatorABCABC(tvec, param, Q)
    )
  }

  K4operatorAABB_factored <- function(tvec, param, B, dvec) {
    r <- .validate_factored_Q(
      tvec, B, dvec, "RSS K4operatorAABB_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    balanced <- .balance_factored_Q(
      B, dvec, "RSS K4operatorAABB_factored"
    )
    B <- balanced$A
    dvec <- balanced$d

    if (r > nrow(B)) {
      Q <- .factor_block_matrix(B, dvec, B)
      base_K4 <- if (summand_K4AABB_safe) {
        summand_K4AABB_factored(tvec, param, B, dvec)
      } else {
        .coordinate_child_K4_AABB(tvec, param, Q)
      }
      return(K4operatorAABB_core(
        tvec,
        param,
        Q,
        base_K4
      ))
    }

    terms <- .factored_contraction_terms(tvec, param, B, dvec)
    dv <- terms$dv

    T3_Q_v <- .T3_AAB_factored(
      summand_cgf, tvec, param, B, dvec, terms$v
    )
    K4_X <- summand_cgf$.private_api$K4operatorAABB_factored(
      tvec, param, B, dvec
    )

    dv$a1 * K4_X +
      dv$a2 * (
        4 * T3_Q_v + terms$trSQ * terms$trSQ + 2 * terms$tr_SQSQ
      ) +
      dv$a3 * (
        2 * terms$trSQ * terms$qmm + 4 * terms$quad_vSv
      ) +
      dv$a4 * (terms$qmm * terms$qmm)
  }

  K3K3operatorAABBCC_factored <- function(tvec, param, B, dvec) {
    r <- .validate_factored_Q(
      tvec, B, dvec, "RSS K3K3operatorAABBCC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    balanced <- .balance_factored_Q(
      B, dvec, "RSS K3K3operatorAABBCC_factored"
    )
    B <- balanced$A
    dvec <- balanced$d

    if (r > nrow(B)) {
      Q <- .factor_block_matrix(B, dvec, B)
      base_T3T3 <- if (summand_K3K3AABBCC_safe) {
        summand_K3K3AABBCC_factored(tvec, param, B, dvec)
      } else {
        .coordinate_child_AABBCC(tvec, param, Q)
      }
      return(K3K3operatorAABBCC_core(
        tvec, param, Q, base_T3T3
      ))
    }

    terms <- .factored_contraction_terms(tvec, param, B, dvec)
    dv <- terms$dv

    h <- terms$trSQ * terms$B_mu + 2 * terms$G_weighted_B_mu
    Q_uU <- as.vector(B %*% (dvec * h))
    T3_Q_v <- .T3_AAB_factored(
      summand_cgf, tvec, param, B, dvec, terms$v
    )
    T3_Q_QuU <- .T3_AAB_factored(
      summand_cgf, tvec, param, B, dvec, Q_uU
    )
    base_T3T3 <-
      summand_cgf$.private_api$K3K3operatorAABBCC_factored(
        tvec, param, B, dvec
      )

    UU <- sum((dvec * h) * h)
    UM <- terms$qmm * sum(h * terms$weighted_B_mu)
    MM <- terms$qmm * terms$qmm * terms$qmm

    (dv$a1 * dv$a1) * base_T3T3 +
      2 * dv$a1 * dv$a2 * T3_Q_QuU +
      2 * dv$a1 * dv$a3 * terms$qmm * T3_Q_v +
      (dv$a2 * dv$a2) * UU +
      2 * dv$a2 * dv$a3 * UM +
      (dv$a3 * dv$a3) * MM
  }

  K3K3operatorABCABC_factored <- function(tvec, param, B, dvec) {
    r <- .validate_factored_Q(
      tvec, B, dvec, "RSS K3K3operatorABCABC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    balanced <- .balance_factored_Q(
      B, dvec, "RSS K3K3operatorABCABC_factored"
    )
    B <- balanced$A
    dvec <- balanced$d

    if (r > nrow(B)) {
      Q <- .factor_block_matrix(B, dvec, B)
      base_T3T3 <- if (summand_K3K3ABCABC_safe) {
        summand_K3K3ABCABC_factored(tvec, param, B, dvec)
      } else {
        .coordinate_child_ABCABC(tvec, param, Q)
      }
      return(K3K3operatorABCABC_core(
        tvec, param, Q, base_T3T3
      ))
    }

    terms <- .factored_contraction_terms(tvec, param, B, dvec)
    dv <- terms$dv

    base_T3T3 <-
      summand_cgf$.private_api$K3K3operatorABCABC_factored(
        tvec, param, B, dvec
      )
    T3_vvv <- summand_cgf$K3operator(
      tvec, param, terms$v, terms$v, terms$v
    )
    T3_M_v <- .T3_QSigmaQ_v_factored(
      summand_cgf, tvec, param, B, dvec, terms$G, terms$v
    )

    sv_Q_sv <- sum(
      (dvec * terms$G_weighted_B_mu) * terms$G_weighted_B_mu
    )
    UU <- 3 * terms$tr_SQSQ * terms$qmm + 6 * sv_Q_sv
    UM <- 3 * terms$quad_vSv * terms$qmm
    MM <- terms$qmm * terms$qmm * terms$qmm

    (dv$a1 * dv$a1) * base_T3T3 +
      6 * dv$a1 * dv$a2 * T3_M_v +
      2 * dv$a1 * dv$a3 * T3_vvv +
      (dv$a2 * dv$a2) * UU +
      2 * dv$a2 * dv$a3 * UM +
      (dv$a3 * dv$a3) * MM
  }
  # Each method selects a factor-space core only while rank does not exceed
  # coordinate dimension, and otherwise uses the dense/coordinate core above.
  # It is therefore safe for a singleton structural wrapper to delegate.
  K4operatorAABB_factored <- .factored_delegate_mark(
    K4operatorAABB_factored, TRUE
  )
  K3K3operatorAABBCC_factored <- .factored_delegate_mark(
    K3K3operatorAABBCC_factored, TRUE
  )
  K3K3operatorABCABC_factored <- .factored_delegate_mark(
    K3K3operatorABCABC_factored, TRUE
  )

  func_T <- function(tvec, param) {
    d <- length(tvec)
    self_object <- get("self", inherits = TRUE)

    Q <- self_object$K2_solve(tvec, param, diag(d))  # Q = K2^{-1}

    K4_AABB <- self_object$K4operatorAABB(tvec, param, Q)
    K3K3_AABBCC <- self_object$K3K3operatorAABBCC(tvec, param, Q)
    K3K3_ABCABC <- self_object$K3K3operatorABCABC(tvec, param, Q)
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
    K2_factor = K2_factor,
    rsim = rsim,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    func_T = func_T,
    ineq_constraint = ineq_constraint,
    op_name = op_name_vec
  )

  contraction_pairs <- list(
    c("K4operatorAABB", "K4operatorAABB_factored"),
    c("K3K3operatorAABBCC", "K3K3operatorAABBCC_factored"),
    c("K3K3operatorABCABC", "K3K3operatorABCABC_factored")
  )
  for (pair in contraction_pairs) {
    if (any(pair %in% extra_names)) {
      cgf_args[[pair[[1L]]]] <- NULL
      cgf_args[[pair[[2L]]]] <- NULL
    }
  }

  contraction_names <- unlist(contraction_pairs, use.names = FALSE)
  protected_names <- setdiff(
    names(cgf_args), c(contraction_names, "func_T")
  )
  conflicting_names <- intersect(extra_names, protected_names)
  if (length(conflicting_names) > 0L) {
    stop(
      "randomlyStoppedSumCGF cannot override generated method(s) through ",
      "'...': ", paste(conflicting_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (any(contraction_names %in% extra_names)) cgf_args$func_T <- NULL

  if ("func_T" %in% extra_names && !is.function(extra_args$func_T)) {
    stop("'func_T' must be a function.", call. = FALSE)
  }

  do.call(createCGF, modifyList(cgf_args, extra_args))
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

  if (is.null(iidReps)) iidReps <- "any"

  # If iidReps="any", block_size is required
  if (identical(iidReps, "any") && is.null(block_size)) {
    stop("randomlyStoppedSumCGF(): iidReps='any' requires a non-NULL 'block_size' ")
  }

  iidReplicatesCGF(cgf = out_cgf, iidReps = iidReps, block_size = block_size)
}
