# -----------------------------------------------------------------------------
# Internal structured K2 factor support
#
# A factor method has signature
#
#   K2_factor(tvec, parameter_vector, A)
#
# and returns additive terms of either list(B = ..., d = ...) or list(S = ...)
# representing
#
#   A K2(t, theta) A' = sum S + sum B diag(d) B'.
#
# The representation lets compositions combine singular covariance
# contributions before requiring the completed covariance to be invertible.
# CGFs without this capability retain their existing dense methods.
# -----------------------------------------------------------------------------

#' @keywords internal
.K2_factor_method <- function(cgf) {
  method <- cgf$additional_methods[["K2_factor"]]
  if (is.function(method)) method else NULL
}

#' @keywords internal
.K2_factor_terminal_method <- function(cgf) {
  factor_method <- .K2_factor_method(cgf)
  use_terminal <- cgf$additional_methods[["K2_factor_terminal"]]
  if (!is.null(factor_method) && is.function(use_terminal) &&
      isTRUE(use_terminal())) {
    factor_method
  } else {
    NULL
  }
}

# A structured solve/logdet pair may be preferred by an enclosing covariance
# update even when a factor representation is also available.  The marker is
# construction-time provenance only.
#' @keywords internal
.K2_structured_pair_mark <- function(method, safe) {
  if (is.null(method)) return(NULL)
  if (!is.function(method)) {
    stop("Internal K2 structured-pair capability requires a function.")
  }
  attr(method, "saddlepoint.K2_structured_pair") <- isTRUE(safe)
  method
}

#' @keywords internal
.K2_structured_pair_is_safe <- function(cgf) {
  capability <- cgf$.private_api$K2_structured_pair
  is.function(capability) && isTRUE(capability())
}

#' @keywords internal
.K2_factor_term <- function(B, d) {
  list(list(B = B, d = as.vector(d)))
}

#' @keywords internal
.K2_dense_term <- function(S) {
  if (is.null(dim(S))) {
    n <- as.integer(sqrt(length(S)))
    if (n * n != length(S)) {
      stop("Dense K2 factor term must be square.", call. = FALSE)
    }
    attr(S, "dim") <- c(n, n)
  }
  list(list(S = S))
}

#' @keywords internal
.K2_factor_scale <- function(terms, scale) {
  lapply(terms, function(term) {
    if (!is.null(term$S)) {
      list(S = term$S * scale)
    } else {
      list(B = term$B, d = term$d * scale)
    }
  })
}

# Build a row-equilibrated Cholesky root of
#
#   S + B diag(d) B'.
#
# sqrt(d) is used only in the numeric forward calculation.  The custom reverse
# rule differentiates B diag(d) B' directly, keeping derivatives defined when
# d is exactly zero or underflows to zero.  The atomic is stateless so an
# invalid optimizer evaluation cannot contaminate a later valid evaluation.
#' @keywords internal
.weighted_gram_chol_atomic <- local({
  forward <- function(x) {
    x <- .rtmb_value_real(x)
    m <- as.integer(x[1L])
    r <- as.integer(x[2L])
    has_dense <- as.integer(x[3L])
    bad <- rep(NaN, m * m + m)

    if (m < 1L || !(has_dense %in% 0:1)) return(bad)

    cursor <- 4L
    S <- base::matrix(0, nrow = m, ncol = m)
    if (has_dense == 1L) {
      S <- base::matrix(x[cursor:(cursor + m * m - 1L)],
                        nrow = m, ncol = m)
      cursor <- cursor + m * m
      if (any(!is.finite(S))) return(bad)
    }

    packed_values <- if (r == 0L) numeric(0) else x[cursor:length(x)]
    packed <- base::matrix(packed_values, nrow = m + 1L, ncol = r)
    d <- base::as.numeric(packed[1L, ])
    B <- base::matrix(
      base::as.numeric(packed[-1L, , drop = FALSE]),
      nrow = m, ncol = r
    )
    if (any(!is.finite(d)) || any(d < 0) || any(!is.finite(B))) {
      return(bad)
    }

    F <- B * rep(sqrt(d), each = m)
    if (any(!is.finite(F))) return(bad)

    # Overflow-safe Euclidean norm of each factor row.
    factor_scale <- if (r == 0L) {
      rep(0, m)
    } else {
      row_max <- apply(abs(F), 1L, max)
      row_divisor <- ifelse(row_max == 0, 1, row_max)
      row_max * sqrt(rowSums((F / row_divisor)^2))
    }

    if (has_dense == 0L) {
      if (r < m || any(!is.finite(factor_scale)) ||
          any(factor_scale <= 0)) return(bad)

      F_scaled <- F / factor_scale
      reference <- base::tcrossprod(F_scaled)
      if (any(!is.finite(reference))) return(bad)
      R_scaled <- tryCatch(base::chol(reference), error = function(e) NULL)

      # Direct-factor QR is a recovery route only when rounding the Gram
      # matrix has hidden positive definiteness.  The ordinary path therefore
      # retains the same inexpensive Gram-plus-Cholesky structure as before.
      if (is.null(R_scaled)) {
        qr_pivoted <- tryCatch(
          base::qr(t(F_scaled), LAPACK = TRUE),
          error = function(e) NULL
        )
        if (is.null(qr_pivoted)) return(bad)
        R_pivoted <- base::qr.R(qr_pivoted, complete = FALSE)
        if (!identical(dim(R_pivoted), c(m, m)) ||
            any(!is.finite(R_pivoted))) return(bad)

        unpivoted_root <- base::matrix(0, nrow = m, ncol = m)
        unpivoted_root[, qr_pivoted$pivot] <- R_pivoted
        qr_canonical <- tryCatch(
          base::qr(unpivoted_root, tol = 0, LAPACK = FALSE),
          error = function(e) NULL
        )
        if (is.null(qr_canonical)) return(bad)
        R_scaled <- base::qr.R(qr_canonical, complete = FALSE)
        signs <- ifelse(base::diag(R_scaled) < 0, -1, 1)
        R_scaled <- R_scaled * rep(signs, times = m)

        # Reject only effective rank loss at machine precision.  This is not
        # a moderate condition-number policy.
        singular_values <- tryCatch(
          base::svd(R_scaled, nu = 0L, nv = 0L)$d,
          error = function(e) NULL
        )
        if (is.null(singular_values) || any(!is.finite(singular_values)) ||
            max(singular_values) <= 0 ||
            min(singular_values) <=
              max(1L, m) * .Machine$double.eps * max(singular_values)) {
          return(bad)
        }
      }
    } else {
      # Validate the completed covariance, not its individual contributions.
      completed_diagonal <- base::diag(S) + factor_scale * factor_scale
      if (any(!is.finite(completed_diagonal)) ||
          any(completed_diagonal <= 0)) return(bad)
      factor_scale <- sqrt(completed_diagonal)

      S_scaled_raw <- base::sweep(S, 1L, factor_scale, "/")
      S_scaled_raw <- base::sweep(S_scaled_raw, 2L, factor_scale, "/")
      if (any(!is.finite(S_scaled_raw))) return(bad)

      # Accept ordinary roundoff skew while rejecting a materially asymmetric
      # object that cannot represent a covariance.  The midpoint is formed
      # without adding two near-overflow entries.
      S_scaled_transpose <- t(S_scaled_raw)
      pair_scale <- pmax(abs(S_scaled_raw), abs(S_scaled_transpose))
      pair_divisor <- pair_scale
      pair_divisor[pair_divisor == 0] <- 1
      S_unit <- S_scaled_raw / pair_divisor
      S_transpose_unit <- S_scaled_transpose / pair_divisor
      S_scaled <- pair_scale * ((S_unit + S_transpose_unit) / 2)
      half_skew <- pair_scale * abs((S_unit - S_transpose_unit) / 2)
      symmetry_tolerance <- sqrt(.Machine$double.eps) / 4
      if (any(!is.finite(half_skew)) ||
          max(rowSums(half_skew)) > symmetry_tolerance) {
        return(bad)
      }

      F_scaled <- if (r == 0L) F else F / factor_scale
      reference <- S_scaled
      if (r > 0L) reference <- reference + base::tcrossprod(F_scaled)
      if (any(!is.finite(reference))) return(bad)

      R_scaled <- tryCatch(base::chol(reference), error = function(e) NULL)
      if (is.null(R_scaled)) return(bad)
    }

    if (any(!is.finite(R_scaled)) || any(base::diag(R_scaled) <= 0) ||
        any(!is.finite(reference))) return(bad)

    # A compact implementation check: the returned root must reconstruct the
    # completed, equilibrated covariance to ordinary double precision.
    reconstruction_error <-
      max(abs(base::crossprod(R_scaled) - reference)) /
      max(abs(reference))
    reconstruction_tolerance <-
      100 * max(1L, m) * .Machine$double.eps
    if (!is.finite(reconstruction_error) ||
        reconstruction_error > reconstruction_tolerance) return(bad)

    c(base::as.vector(R_scaled), factor_scale)
  }

  reverse <- function(x, out, out_bar) {
    x <- RTMB::AD(x)
    out <- RTMB::AD(out)
    out_bar <- RTMB::AD(out_bar)

    m <- as.integer((sqrt(1 + 4 * length(out)) - 1) / 2)
    payload_length <- length(x) - 3L
    has_dense <- as.integer(payload_length %% (m + 1L) == 1L)
    r <- as.integer((payload_length - has_dense * m * m) / (m + 1L))
    R_index <- seq_len(m * m)
    scale_index <- m * m + seq_len(m)
    R <- matrix(out[R_index], nrow = m, ncol = m)
    R_bar <- matrix(out_bar[R_index], nrow = m, ncol = m)
    scale <- as.vector(out[scale_index])
    scale_bar <- as.vector(out_bar[scale_index])

    cursor <- 4L
    if (has_dense == 1L) {
      S_index <- cursor:(cursor + m * m - 1L)
      cursor <- cursor + m * m
    } else {
      S_index <- integer(0)
    }
    packed_values <- if (r == 0L) x[integer(0)] else x[cursor:length(x)]
    packed <- matrix(packed_values, nrow = m + 1L, ncol = r)
    d <- as.vector(packed[1L, ])
    B <- packed[-1L, , drop = FALSE]

    # Reverse the canonical Cholesky of E, then reverse K = D E D.
    L <- t(R)
    M <- t(L) %*% t(R_bar)
    lower_mask <- lower.tri(matrix(0, m, m), diag = FALSE)
    P <- M * lower_mask + 0.5 * M * diag(m)
    G <- 0.5 * (P + t(P))
    E_bar <- solve(t(L), G %*% solve(L, diag(m)))
    E_bar <- 0.5 * (E_bar + t(E_bar))
    E <- crossprod(R)

    E_weighted <- E_bar * E
    attr(E_weighted, "dim") <- c(m, m)
    scale_bar <- scale_bar -
      (
        as.vector(E_weighted %*% rep(1, m)) +
          as.vector(t(E_weighted) %*% rep(1, m))
      ) / scale
    K_bar <- E_bar / rep(scale, times = m)
    K_bar <- K_bar / rep(scale, each = m)
    diagonal_adjustment <- matrix(0, nrow = m, ncol = m) * scale_bar[1L]
    diagonal_adjustment[cbind(seq_len(m), seq_len(m))] <-
      scale_bar / (2 * scale)
    K_bar <- K_bar + diagonal_adjustment
    K_bar <- 0.5 * (K_bar + t(K_bar))

    packed_bar <- if (r > 0L) {
      KB <- K_bar %*% B
      B_bar <- 2 * KB * rep(d, each = m)
      d_bar <- colSums(B * KB)
      as.vector(rbind(d_bar, B_bar))
    } else {
      x[integer(0)]
    }

    prefix_bar <- c(0 * x[1L], 0 * x[2L], 0 * x[3L])
    dense_bar <- if (has_dense == 1L) {
      as.vector(K_bar + 0 * x[S_index])
    } else {
      NULL
    }
    c(prefix_bar, dense_bar, packed_bar)
  }

  RTMB::ADjoint(forward, reverse, name = "weighted_gram_chol")
})

#' @keywords internal
.weighted_gram_chol <- local({
  atomic <- .weighted_gram_chol_atomic
  function(x, S = NULL) {
    m <- nrow(x) - 1L
    r <- ncol(x)
    x_vec <- as.vector(x)
    S_vec <- if (is.null(S)) numeric(0) else as.vector(S)
    payload <- RTMB::AD(numeric(length(x_vec) + length(S_vec) + 3L))
    payload[1:3] <- c(m, r, !is.null(S))
    if (length(S_vec) > 0L) {
      payload[3L + seq_along(S_vec)] <- S_vec
    }
    payload[3L + length(S_vec) + seq_along(x_vec)] <- x_vec
    out <- atomic(payload)
    R <- out[seq_len(m * m)]
    attr(R, "dim") <- c(m, m)
    list(R = R, scale = out[m * m + seq_len(m)])
  }
})

# Factor a dense SPD matrix through the same stateless terminal used by
# structured K2 factors.  This is used when a valid thin square root must be
# propagated to higher-order contractions.
#' @keywords internal
.K2_dense_spd_factor <- function(Q, normalize = TRUE) {
  Q_dim <- dim(Q)
  if (is.null(Q_dim) && length(Q) == 1L) {
    attr(Q, "dim") <- c(1L, 1L)
    Q_dim <- c(1L, 1L)
  }
  if (length(Q_dim) != 2L || Q_dim[1L] < 1L || Q_dim[1L] != Q_dim[2L]) {
    stop("Dense K2 factor input must be a non-empty square matrix.",
         call. = FALSE)
  }

  n <- Q_dim[1L]
  packed <- matrix(numeric(0), nrow = n + 1L, ncol = 0L)
  factorization <- .weighted_gram_chol(packed, Q)
  if (!inherits(factorization$R, "advector") &&
      any(!is.finite(c(factorization$R, factorization$scale)))) {
    stop("Matrix is singular, indefinite, or numerically invalid.",
         call. = FALSE)
  }

  B <- factorization$scale * t(factorization$R)
  if (!normalize) return(list(B = B, d = rep(1, n)))

  diagonal <- factorization$scale * diag(factorization$R)
  list(
    B = B %*% diag(1 / diagonal, nrow = n, ncol = n),
    d = diagonal * diagonal
  )
}

# Check only the componentwise backward residual of the completed scaled
# solve.  This is an implementation sanity check, not a guarantee of forward
# or derivative accuracy.  Its pullback is the identity on successful values.
#' @keywords internal
.K2_residual_check_atomic <- local({
  forward <- function(x) {
    x <- .rtmb_value_real(x)
    m <- as.integer(x[1L])
    nrhs <- as.integer(x[2L])
    solution_length <- m * nrhs
    bad <- rep(NaN, solution_length + 1L)
    expected_length <- 2L + m * m + 3L * solution_length
    if (m < 1L || nrhs < 1L || length(x) != expected_length ||
        any(!is.finite(x))) return(bad)

    cursor <- 3L
    solution <- base::matrix(
      x[cursor:(cursor + solution_length - 1L)],
      nrow = m, ncol = nrhs
    )
    cursor <- cursor + solution_length
    R <- base::matrix(x[cursor:(cursor + m * m - 1L)],
                      nrow = m, ncol = m)
    cursor <- cursor + m * m
    w <- base::matrix(x[cursor:(cursor + solution_length - 1L)],
                      nrow = m, ncol = nrhs)
    cursor <- cursor + solution_length
    z <- base::matrix(x[cursor:(cursor + solution_length - 1L)],
                      nrow = m, ncol = nrhs)
    if (any(base::diag(R) <= 0)) return(bad)

    Rw <- R %*% w
    residual <- z - t(R) %*% Rw
    # |R'R| |w| is bounded by |R'| |R| |w|.  Using the triangular
    # factors avoids rebuilding the dense covariance for every right-hand
    # side while retaining a standard componentwise backward-error scale.
    denominator <- abs(z) + t(abs(R)) %*% (abs(R) %*% abs(w))
    active <- denominator > 0
    if (any(!active & residual != 0)) return(bad)
    backward_error <- if (any(active)) {
      max(abs(residual[active]) / denominator[active])
    } else {
      0
    }
    backward_tolerance <- 64 * max(1L, m) * .Machine$double.eps
    if (!is.finite(backward_error) ||
        backward_error > backward_tolerance) return(bad)

    c(as.vector(solution), 1)
  }

  reverse <- function(x, out, out_bar) {
    x <- RTMB::AD(x)
    out <- RTMB::AD(out)
    out_bar <- RTMB::AD(out_bar)
    solution_length <- length(out) - 1L
    solution_index <- 2L + seq_len(solution_length)
    status <- out[solution_length + 1L]
    status_derivative <- status / status - 1
    solution_bar <- out_bar[seq_len(solution_length)] * status +
      out_bar[solution_length + 1L] * status_derivative
    c(
      0 * x[1L],
      0 * x[2L],
      solution_bar + 0 * x[solution_index],
      0 * x[-c(1L, 2L, solution_index)]
    )
  }

  RTMB::ADjoint(forward, reverse, name = "K2_residual_check")
})

#' @keywords internal
.K2_residual_check <- local({
  atomic <- .K2_residual_check_atomic
  function(solution, R, w, z) {
    m <- nrow(R)
    nrhs <- if (is.null(dim(solution))) 1L else ncol(solution)
    solution_vec <- as.vector(solution)
    payload <- RTMB::AD(numeric(
      2L + m * m + 3L * length(solution_vec)
    ))
    payload[1:2] <- c(m, nrhs)
    cursor <- 3L
    payload[cursor:(cursor + length(solution_vec) - 1L)] <- solution_vec
    cursor <- cursor + length(solution_vec)
    payload[cursor:(cursor + m * m - 1L)] <- as.vector(R)
    cursor <- cursor + m * m
    payload[cursor:(cursor + length(solution_vec) - 1L)] <- as.vector(w)
    cursor <- cursor + length(solution_vec)
    payload[cursor:(cursor + length(solution_vec) - 1L)] <- as.vector(z)
    atomic_out <- atomic(payload)
    out <- atomic_out[seq_along(solution_vec)]
    if (!is.null(dim(solution))) attr(out, "dim") <- dim(solution)
    out
  }
})

#' @keywords internal
.K2_factor_chol <- function(K2_factor, tvec, parameter_vector) {
  n <- length(tvec)
  terms <- K2_factor(tvec, parameter_vector, diag(n))
  if (!is.list(terms) || length(terms) == 0L) {
    stop("K2 factorization supplied no covariance terms.", call. = FALSE)
  }

  dense_terms <- list()
  packed_terms <- list()
  for (term in terms) {
    if (!is.null(term$S)) {
      if (!identical(dim(term$S), c(n, n))) {
        stop("Invalid K2 factor term dimensions.", call. = FALSE)
      }
      dense_terms[[length(dense_terms) + 1L]] <- term$S
    } else if (!is.null(term$B) && !is.null(term$d)) {
      di <- as.vector(term$d)
      if (nrow(term$B) != n || ncol(term$B) != length(di)) {
        stop("Invalid K2 factor term dimensions.", call. = FALSE)
      }
      Bi <- term$B
      if (inherits(Bi, "denseMatrix") && !inherits(Bi, "adsparse")) {
        Bi <- as.matrix(Bi)
      }
      packed_terms[[length(packed_terms) + 1L]] <- rbind(di, Bi)
    } else {
      stop("Invalid K2 factor term.", call. = FALSE)
    }
  }

  S <- if (length(dense_terms) == 0L) {
    NULL
  } else if (length(dense_terms) == 1L) {
    dense_terms[[1L]]
  } else {
    Reduce(`+`, dense_terms)
  }
  packed <- if (length(packed_terms) == 0L) {
    matrix(numeric(0), nrow = n + 1L, ncol = 0L)
  } else if (length(packed_terms) == 1L) {
    packed_terms[[1L]]
  } else {
    do.call(cbind, packed_terms)
  }
  factorization <- .weighted_gram_chol(packed, S)

  if (!inherits(factorization$R, "advector") &&
      any(!is.finite(c(factorization$R, factorization$scale)))) {
    stop(
      "K2 is singular, indefinite, or numerically invalid at this evaluation.",
      call. = FALSE
    )
  }
  factorization
}

#' @keywords internal
.K2_factor_solve <- function(K2_factor, tvec, parameter_vector, rhs) {
  factorization <- .K2_factor_chol(K2_factor, tvec, parameter_vector)
  R <- factorization$R
  z <- rhs / factorization$scale
  w <- solve(R, solve(t(R), z))
  solution <- w / factorization$scale
  solution <- .K2_residual_check(solution, R, w, z)
  if (!inherits(solution, "advector") && any(!is.finite(solution))) {
    stop("K2 solve failed its numerical residual check.", call. = FALSE)
  }
  solution
}

#' @keywords internal
.K2_factor_logdet <- function(K2_factor, tvec, parameter_vector) {
  factorization <- .K2_factor_chol(K2_factor, tvec, parameter_vector)
  R <- factorization$R
  value <- 2 * sum(log(diag(R))) +
    2 * sum(log(factorization$scale))
  if (!inherits(value, "advector") && any(!is.finite(value))) {
    stop("K2 log-determinant is not finite.", call. = FALSE)
  }
  value
}
