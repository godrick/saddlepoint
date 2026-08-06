# -----------------------------------------------------------------------------
# Internal K2 factor support
#
# A certified factor method has signature
#
#   K2_factor(tvec, parameter_vector, A)
#
# and returns additive terms of either list(B = ..., d = ...) or list(S = ...)
# representing
#
#   A K2(t, theta) A' = sum S + sum B diag(d) B'.
#
# This is deliberately private to the numerical K2 solve/logdet path.  CGFs that
# do not supply the capability continue to use the existing dense methods.
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
  if (!is.null(factor_method) && is.function(use_terminal) && isTRUE(use_terminal())) {
    factor_method
  } else {
    NULL
  }
}

# A structured solve/logdet pair may be preferred by an enclosing covariance
# update even when a factor representation is also available.  The marker is
# construction-time provenance only: package wrappers propagate it
# mechanically, while an untagged pair supplied directly to createCGF() is
# authoritative by default.
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

# The forward factorization sees sqrt(d), but sqrt(d) is intentionally not part
# of the AD tape.  The reverse rule differentiates B diag(d) B' directly, so
# derivatives remain defined at d == 0 (including underflow to zero).  Each
# evaluation is stateless and initializes a fresh invalid result.
#' @keywords internal
.weighted_gram_chol_atomic <- local({
  forward <- function(x) {
    x <- .rtmb_value_real(x)
    m <- as.integer(x[1L])
    r <- as.integer(x[2L])
    has_dense <- as.integer(x[3L])
    bad <- rep(NaN, m * m + 2L * m)

    if (m < 1L || !(has_dense %in% 0:1)) return(bad)

    cursor <- 4L
    S <- base::matrix(0, nrow = m, ncol = m)
    if (has_dense == 1L) {
      S <- base::matrix(x[cursor:(cursor + m * m - 1L)], nrow = m, ncol = m)
      cursor <- cursor + m * m
      if (any(!is.finite(S))) return(bad)
    }

    packed_values <- if (r == 0L) numeric(0) else x[cursor:length(x)]
    packed <- base::matrix(packed_values, nrow = m + 1L, ncol = r)
    d <- base::as.numeric(packed[1L, ])
    B <- base::matrix(
      base::as.numeric(packed[-1L, , drop = FALSE]),
      nrow = m,
      ncol = r
    )
    if (any(!is.finite(d)) || any(d < 0) || any(!is.finite(B))) return(bad)

    F <- B * rep(sqrt(d), each = m)
    if (any(!is.finite(F))) return(bad)

    # Overflow-safe row Euclidean norms.
    factor_scale <- if (r == 0L) {
      rep(0, m)
    } else {
      row_max <- apply(abs(F), 1L, max)
      row_divisor <- ifelse(row_max == 0, 1, row_max)
      row_max * sqrt(rowSums((F / row_divisor)^2))
    }

    if (has_dense == 0L) {
      # Form the Gram matrix only after row equilibration.  Certify its
      # dot-product roundoff with a rowwise a-priori bound.  This needs only
      # O(m*r) work, unlike forming a second m-by-m Gram product, and it also
      # detects the shared rounding loss that two algebraically equivalent
      # BLAS products can reproduce identically.
      if (r < m || any(!is.finite(factor_scale)) || any(factor_scale <= 0)) {
        return(bad)
      }
      F_scaled <- F / factor_scale
      factor_reference <- base::tcrossprod(F_scaled)
      if (any(!is.finite(factor_reference))) return(bad)
      reference <- 0.5 * (factor_reference + t(factor_reference))
      active_columns <- colSums(abs(F_scaled)) > 0
      effective_width <- sum(active_columns)
      gram_operations <- effective_width + 8L
      gram_gamma <- gram_operations * .Machine$double.eps /
        (1 - gram_operations * .Machine$double.eps)
      if (!is.finite(gram_gamma) || gram_gamma >= 1) return(bad)
      formation_discrepancy <- if (effective_width == 0L) {
        rep(0, m)
      } else {
        active_F <- abs(F_scaled[, active_columns, drop = FALSE])
        gram_gamma * as.vector(active_F %*% colSums(active_F))
      }

      # Cholesky is the ordinary route.  Direct-factor QR is only a recovery
      # path when the rounded Gram matrix has lost positive definiteness.  The
      # stricter covariance-level condition and accuracy checks belong to the
      # solve/logdet terminal below; this root builder is also used by mapped
      # higher-order contractions, which do not require an inverse.
      R_scaled <- tryCatch(base::chol(reference), error = function(e) NULL)
      if (is.null(R_scaled)) {
        qr_pivoted <- tryCatch(
          base::qr(t(F_scaled), LAPACK = TRUE),
          error = function(e) NULL
        )
        if (is.null(qr_pivoted)) return(bad)
        R_pivoted <- base::qr.R(qr_pivoted, complete = FALSE)
        if (!identical(dim(R_pivoted), c(m, m)) ||
            any(!is.finite(R_pivoted))) return(bad)

        # Undo LAPACK's coordinate pivot before returning the ordinary-
        # coordinate positive-diagonal root.  A discrepancy from the rounded
        # Gram matrix is retained below and certified by every solve.
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
      # Validate a mixed contribution only after completing it.  A dense term
      # can have a roundoff-sized negative eigenvalue while the final covariance
      # is safely SPD; neither the term nor the completed sum is projected.
      S_scale <- max(abs(S))
      symmetry_tolerance <-
        10 * max(1L, m) * .Machine$double.eps * S_scale
      if (!is.finite(S_scale) ||
          max(abs(S - t(S))) > symmetry_tolerance) return(bad)
      S <- 0.5 * (S + t(S))

      completed_diagonal <- base::diag(S) + factor_scale * factor_scale
      if (any(!is.finite(completed_diagonal)) ||
          any(completed_diagonal <= 0)) return(bad)
      factor_scale <- sqrt(completed_diagonal)

      S_scaled <- base::sweep(S, 1L, factor_scale, "/")
      S_scaled <- base::sweep(S_scaled, 2L, factor_scale, "/")
      F_scaled <- if (r == 0L) F else F / factor_scale
      factor_reference <- S_scaled
      if (r > 0L) {
        factor_reference <- factor_reference + base::tcrossprod(F_scaled)
      }
      if (any(!is.finite(factor_reference))) return(bad)
      reference <- 0.5 * (factor_reference + t(factor_reference))
      active_columns <- if (r == 0L) {
        logical(0)
      } else {
        colSums(abs(F_scaled)) > 0
      }
      effective_width <- sum(active_columns)
      gram_operations <- effective_width + 8L
      gram_gamma <- gram_operations * .Machine$double.eps /
        (1 - gram_operations * .Machine$double.eps)
      if (!is.finite(gram_gamma) || gram_gamma >= 1) return(bad)
      factor_roundoff <- if (effective_width == 0L) {
        rep(0, m)
      } else {
        active_F <- abs(F_scaled[, active_columns, drop = FALSE])
        gram_gamma * as.vector(active_F %*% colSums(active_F))
      }
      # Scaling S and adding it to the factor Gram each contribute only
      # m-by-m rounding.  Include both in the same rowwise certificate.
      dense_gamma <- (2L * max(1L, m) + 4L) * .Machine$double.eps
      formation_discrepancy <- factor_roundoff +
        dense_gamma * rowSums(abs(S_scaled)) +
        .Machine$double.eps * rowSums(abs(reference))

      R_scaled <- tryCatch(base::chol(reference), error = function(e) NULL)
      if (is.null(R_scaled) || any(!is.finite(R_scaled))) return(bad)
    }

    if (any(!is.finite(R_scaled)) || any(base::diag(R_scaled) <= 0)) return(bad)
    reconstruction_delta <- abs(base::crossprod(R_scaled) - reference)
    root_operations <- max(1L, m) + 8L
    root_gamma <- root_operations * .Machine$double.eps /
      (1 - root_operations * .Machine$double.eps)
    if (!is.finite(root_gamma) || root_gamma >= 1) return(bad)
    abs_root <- abs(R_scaled)
    root_roundoff <- root_gamma * as.vector(
      t(abs_root) %*% rowSums(abs_root)
    )
    if (any(!is.finite(root_roundoff))) return(bad)
    reconstruction_discrepancy <-
      rowSums(reconstruction_delta) + root_roundoff
    reconstruction_error <- max(reconstruction_delta) / max(abs(reference))
    reconstruction_tolerance <-
      100 * max(1L, m) * .Machine$double.eps
    if (!is.finite(reconstruction_error) ||
        reconstruction_error > reconstruction_tolerance) return(bad)

    # Carry rowwise factor-formation and root-reconstruction loss into every
    # solve certificate.  Keeping rows separate prevents one local discrepancy
    # from being charged to every coordinate.
    solve_discrepancy <- formation_discrepancy + reconstruction_discrepancy
    c(
      base::as.vector(R_scaled), factor_scale,
      solve_discrepancy
    )
  }

  reverse <- function(x, out, out_bar) {
    x <- RTMB::AD(x)
    out <- RTMB::AD(out)
    out_bar <- RTMB::AD(out_bar)

    m <- as.integer(sqrt(length(out) + 1L) - 1L)
    payload_length <- length(x) - 3L
    has_dense <- as.integer(payload_length %% (m + 1L) == 1L)
    r <- as.integer((payload_length - has_dense * m * m) / (m + 1L))
    R_index <- seq_len(m * m)
    scale_index <- m * m + seq_len(m)
    discrepancy_index <- m * m + m + seq_len(m)
    R <- matrix(out[R_index], nrow = m, ncol = m)
    R_bar <- matrix(out_bar[R_index], nrow = m, ncol = m)
    scale <- as.vector(out[scale_index])
    scale_bar <- as.vector(out_bar[scale_index])
    discrepancy_bar <- out_bar[discrepancy_index]

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

    # Reverse the equilibrated Cholesky, then K = D E D with
    # D=sqrt(diag(K)).  All triangular solves use the equilibrated root.
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

    prefix_bar <- c(0 * x[1L], 0 * x[2L], 0 * x[3L]) +
      0 * sum(discrepancy_bar)
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
    list(
      R = R,
      scale = out[m * m + seq_len(m)],
      covariance_discrepancy = out[m * m + m + seq_len(m)]
    )
  }
})

# Factor a dense SPD matrix through the same stateless, certified terminal used
# by structured K2 factors.  This is intentionally private and is used only
# where a valid thin square root must be propagated to higher-order operators.
#' @keywords internal
.K2_dense_spd_factor <- function(Q, normalize = TRUE) {
  Q_dim <- dim(Q)
  if (is.null(Q_dim) && length(Q) == 1L) {
    attr(Q, "dim") <- c(1L, 1L)
    Q_dim <- c(1L, 1L)
  }
  if (length(Q_dim) != 2L || Q_dim[1L] < 1L || Q_dim[1L] != Q_dim[2L]) {
    stop("Dense K2 factor input must be a non-empty square matrix.", call. = FALSE)
  }

  n <- Q_dim[1L]
  packed <- matrix(numeric(0), nrow = n + 1L, ncol = 0L)
  factorization <- .weighted_gram_chol(packed, Q)
  if (!inherits(factorization$R, "advector") &&
      any(!is.finite(c(factorization$R, factorization$scale)))) {
    stop(
      "Matrix is singular, indefinite, or numerically ill-conditioned.",
      call. = FALSE
    )
  }

  B <- factorization$scale * t(factorization$R)
  if (!normalize) return(list(B = B, d = rep(1, n)))

  diagonal <- factorization$scale * diag(factorization$R)
  list(
    B = B %*% diag(1 / diagonal, nrow = n, ncol = n),
    d = diagonal * diagonal
  )
}

# A root can be useful for a higher-order contraction even when its covariance
# is too ill-conditioned to invert accurately.  Keep the stronger condition
# and factor-formation checks in this terminal-only guard.  It returns the
# condition information used by the solve certificate, avoiding another SVD.
#' @keywords internal
.K2_terminal_root_guard_atomic <- local({
  forward <- function(x) {
    x <- .rtmb_value_real(x)
    derivative_mode <- as.integer(x[length(x)])
    payload <- x[-length(x)]
    m <- as.integer((sqrt(1 + 4 * length(payload)) - 1) / 2)
    bad <- rep(NaN, m * m + 4L)
    if (!(derivative_mode %in% 0:1) || m < 1L ||
        m * (m + 1L) != length(payload) ||
        any(!is.finite(x))) return(bad)
    R <- base::matrix(payload[seq_len(m * m)], nrow = m, ncol = m)
    covariance_discrepancy <- payload[m * m + seq_len(m)]
    if (any(covariance_discrepancy < 0)) return(bad)

    accuracy_floor <- 4 * max(1L, m) * sqrt(.Machine$double.eps)
    derivative_condition_floor <- if (derivative_mode == 1L) {
      .Machine$double.eps^(1 / 4)
    } else {
      0
    }
    condition_floor <- max(accuracy_floor, derivative_condition_floor)
    discrepancy_tolerance <- sqrt(.Machine$double.eps) / 4

    # Most equilibrated roots are comfortably conditioned.  Before paying for
    # a full SVD, try a conservative spectral proof based only on positive
    # products.  For E = R'R with eigenvalues lambda_1 <= ... <= lambda_m,
    #
    #   lambda_1 >= det(E) / (trace(E) / (m - 1))^(m - 1),
    #   lambda_m <= trace(E).
    #
    # Every floating-point quantity below is rounded outwards.  Relative
    # product bounds are used only in the normal range; underflow, overflow,
    # or an inconclusive bound falls through to the unchanged SVD route.
    gamma_n <- function(k) {
      ku <- k * .Machine$double.eps
      if (!is.finite(ku) || ku >= 1) Inf else ku / (1 - ku)
    }
    tiny <- .Machine$double.xmin
    R_values <- base::as.vector(R)
    squared_values <- R_values * R_values
    diagonal_squared <- base::diag(R) * base::diag(R)
    fast_spectral_proof <- FALSE
    lambda_lower <- 0
    rho_lower <- 0

    if (all(is.finite(squared_values)) &&
        all(is.finite(diagonal_squared)) &&
        all(diagonal_squared >= tiny)) {
      trace_gamma <- gamma_n(length(squared_values) + 16L)
      trace_upper <- (
        sum(squared_values) + length(squared_values) * tiny
      ) / (1 - trace_gamma)

      determinant_path <- cumprod(diagonal_squared)
      determinant_observed <- determinant_path[length(determinant_path)]
      determinant_gamma <- gamma_n(4L * m + 16L)
      determinant_lower <- if (
        all(is.finite(determinant_path)) &&
          all(determinant_path >= tiny) &&
          is.finite(determinant_gamma)
      ) {
        determinant_observed * (1 - determinant_gamma)
      } else {
        0
      }

      if (is.finite(trace_upper) && trace_upper > 0 &&
          determinant_lower > 0) {
        if (m == 1L) {
          lambda_lower <- determinant_lower
        } else {
          mean_upper <- (trace_upper / (m - 1L)) /
            (1 - 8 * .Machine$double.eps)
          denominator_observed <- prod(rep(mean_upper, m - 1L))
          power_gamma <- gamma_n(m + 16L)
          denominator_upper <- if (
            is.finite(denominator_observed) &&
              denominator_observed >= tiny && is.finite(power_gamma)
          ) {
            denominator_observed / (1 - power_gamma)
          } else {
            Inf
          }
          quotient <- determinant_lower / denominator_upper
          lambda_lower <- if (is.finite(quotient) && quotient > 0) {
            quotient * (1 - 8 * .Machine$double.eps)
          } else {
            0
          }
        }

        quotient <- lambda_lower / trace_upper
        rho_lower <- if (is.finite(quotient) && quotient > 0) {
          quotient * (1 - 8 * .Machine$double.eps)
        } else {
          0
        }
        quotient <- max(covariance_discrepancy) / lambda_lower
        discrepancy_ratio_upper <- if (
          is.finite(quotient) && lambda_lower > 0
        ) {
          quotient / (1 - 8 * .Machine$double.eps) + 8 * tiny
        } else {
          Inf
        }
        fast_spectral_proof <-
          is.finite(rho_lower) && rho_lower > condition_floor &&
          is.finite(discrepancy_ratio_upper) &&
          discrepancy_ratio_upper < 1 &&
          discrepancy_ratio_upper <= discrepancy_tolerance
      }
    }

    if (fast_spectral_proof) {
      return(c(
        base::as.vector(R), rho_lower, lambda_lower,
        1, 1
      ))
    }

    singular_values <- tryCatch(
      base::svd(R, nu = 0L, nv = 0L)$d,
      error = function(e) NULL
    )
    if (is.null(singular_values) || any(!is.finite(singular_values)) ||
        max(singular_values) <= 0) return(bad)
    rho_K <- (min(singular_values) / max(singular_values))^2
    if (!is.finite(rho_K) || rho_K <= accuracy_floor) return(bad)

    # Second derivatives of a solve contain up to three inverse factors.  When
    # this terminal is recorded on an AD tape, reject reciprocal conditions
    # below u^(1/4), where the corresponding generic cubic roundoff
    # amplification can exceed u^(1/4) (about 1e-4 in double precision).
    # Numeric or parameter-independent covariance terms retain the ordinary
    # condition rule above.
    if (rho_K <= derivative_condition_floor) {
      return(bad)
    }
    lambda_min <- min(singular_values)^2
    discrepancy_ratio <- max(covariance_discrepancy) / lambda_min
    if (!is.finite(lambda_min) || lambda_min <= 0 ||
        !is.finite(discrepancy_ratio) || discrepancy_ratio >= 1 ||
        discrepancy_ratio > discrepancy_tolerance) return(bad)
    c(base::as.vector(R), rho_K, lambda_min, 0, 1)
  }

  reverse <- function(x, out, out_bar) {
    x <- RTMB::AD(x)
    out <- RTMB::AD(out)
    out_bar <- RTMB::AD(out_bar)
    m <- as.integer(sqrt(length(out) - 4L))
    status_index <- m * m + 4L
    status <- out[status_index]
    status_derivative <- status / status - 1
    root_bar <- out_bar[seq_len(m * m)] * status +
      out_bar[status_index] * status_derivative
    c(
      root_bar + 0 * x[seq_len(m * m)],
      0 * x[(m * m + 1L):length(x)]
    )
  }

  RTMB::ADjoint(forward, reverse, name = "K2_terminal_root_guard")
})

#' @keywords internal
.K2_terminal_root_guard <- local({
  atomic <- .K2_terminal_root_guard_atomic
  function(R, covariance_discrepancy, covariance_is_ad) {
    m <- nrow(R)
    derivative_mode <- as.integer(isTRUE(covariance_is_ad))
    out <- atomic(c(
      as.vector(R), covariance_discrepancy, derivative_mode
    ))
    guarded_R <- out[seq_len(m * m)]
    attr(guarded_R, "dim") <- dim(R)
    list(
      R = guarded_R,
      rho = out[m * m + 1L],
      lambda_min = out[m * m + 2L],
      condition_is_bound = out[m * m + 3L]
    )
  }
})

# Certify the actual scaled solve at each primal evaluation.  The atomic is the
# identity on every valid evaluation, so its pullback is exactly the identity
# and adds no derivative approximation or mutable cache.  As with any local
# solve, this cannot inspect cancellation introduced later by an enclosing AD
# tangent or adjoint contraction.
#' @keywords internal
.K2_solution_guard_atomic <- local({
  exact_root_condition <- function(R, covariance_discrepancy,
                                   condition_tolerance) {
    singular_values <- tryCatch(
      base::svd(R, nu = 0L, nv = 0L)$d,
      error = function(e) NULL
    )
    if (is.null(singular_values) || any(!is.finite(singular_values)) ||
        max(singular_values) <= 0) return(NULL)

    rho_K <- (min(singular_values) / max(singular_values))^2
    lambda_min <- min(singular_values)^2
    discrepancy_ratio <- max(covariance_discrepancy) / lambda_min
    discrepancy_tolerance <- sqrt(.Machine$double.eps) / 4
    if (!is.finite(rho_K) || rho_K <= condition_tolerance ||
        !is.finite(lambda_min) || lambda_min <= 0 ||
        !is.finite(discrepancy_ratio) || discrepancy_ratio >= 1 ||
        discrepancy_ratio > discrepancy_tolerance) return(NULL)

    list(rho = rho_K, lambda_min = lambda_min)
  }

  forward <- function(x) {
    x <- .rtmb_value_real(x)
    m <- as.integer(x[1L])
    nrhs <- as.integer(x[2L])
    solution_length <- m * nrhs
    bad <- rep(NaN, solution_length + 1L)
    expected_length <- 5L + m * m + 3L * solution_length + 2L * m
    if (m < 1L || nrhs < 1L || length(x) != expected_length ||
        any(!is.finite(x))) return(bad)

    cursor <- 3L
    solution <- base::matrix(
      x[cursor:(cursor + solution_length - 1L)],
      nrow = m,
      ncol = nrhs
    )
    cursor <- cursor + solution_length
    R <- base::matrix(x[cursor:(cursor + m * m - 1L)], nrow = m, ncol = m)
    cursor <- cursor + m * m
    w <- base::matrix(x[cursor:(cursor + solution_length - 1L)],
                      nrow = m, ncol = nrhs)
    cursor <- cursor + solution_length
    z <- base::matrix(x[cursor:(cursor + solution_length - 1L)],
                      nrow = m, ncol = nrhs)
    cursor <- cursor + solution_length
    scale <- x[cursor:(cursor + m - 1L)]
    cursor <- cursor + m
    covariance_discrepancy <- x[cursor:(cursor + m - 1L)]
    cursor <- cursor + m
    rho_K <- x[cursor]
    lambda_min <- x[cursor + 1L]
    condition_is_bound <- x[cursor + 2L]
    if (any(scale <= 0) || any(covariance_discrepancy < 0)) return(bad)
    if (!(condition_is_bound %in% 0:1)) return(bad)
    condition_tolerance <- 4 * max(1L, m) * sqrt(.Machine$double.eps)
    if (!is.finite(rho_K) || rho_K <= condition_tolerance ||
        !is.finite(lambda_min) || lambda_min <= 0) return(bad)

    E <- base::crossprod(R)
    E_norm <- max(rowSums(abs(E)))
    residual <- z - E %*% w
    backward_tolerance <- 64 * max(1L, m) * .Machine$double.eps
    roundoff_floor <- 8 * max(1L, m) * .Machine$double.eps
    columns_to_certify <- integer(0)
    eta <- numeric(nrhs)
    for (j in seq_len(nrhs)) {
      denominator <- E_norm * max(abs(w[, j])) + max(abs(z[, j]))
      eta[j] <- if (denominator == 0) {
        0
      } else {
        max(abs(residual[, j])) / denominator
      }
      if (!is.finite(eta[j]) || eta[j] > backward_tolerance) return(bad)

      if (all(z[, j] == 0)) {
        if (any(solution[, j] != 0)) return(bad)
      } else if (max(abs(solution[, j])) == 0) {
        return(bad)
      } else {
        columns_to_certify <- c(columns_to_certify, j)
      }
    }

    relative_error_proxy <- pmax(eta, roundoff_floor) / rho_K
    if (any(!is.finite(relative_error_proxy)) ||
        any(relative_error_proxy >= 1)) {
      if (condition_is_bound != 1) return(bad)
      exact <- exact_root_condition(
        R, covariance_discrepancy, condition_tolerance
      )
      if (is.null(exact)) return(bad)
      rho_K <- exact$rho
      lambda_min <- exact$lambda_min
      condition_is_bound <- 0
      relative_error_proxy <- pmax(eta, roundoff_floor) / rho_K
      if (any(!is.finite(relative_error_proxy)) ||
          any(relative_error_proxy >= 1)) return(bad)
    }

    # First try the inexpensive triangular roundoff bound.  It is deliberately
    # conservative, but when it succeeds no additional inverse is needed.  If
    # it is inconclusive, use a componentwise a-posteriori residual bound below
    # rather than rejecting an accurately solved ordinary system.
    if (length(columns_to_certify) > 0L) {
      comparison_gamma <- 2 * max(1L, m) * .Machine$double.eps
      if (!is.finite(comparison_gamma) || comparison_gamma >= 1) return(bad)

      w_checked <- w[, columns_to_certify, drop = FALSE]
      z_checked <- z[, columns_to_certify, drop = FALSE]
      solution_checked <- solution[, columns_to_certify, drop = FALSE]
      abs_R <- abs(R)
      Rw <- R %*% w_checked
      raw_residual <- z_checked - t(R) %*% Rw
      residual_bound <- abs(raw_residual) +
        covariance_discrepancy %o%
          apply(abs(w_checked), 2L, max)
      if (any(!is.finite(residual_bound))) return(bad)

      first_magnitude <- base::matrix(
        0, nrow = m, ncol = length(columns_to_certify)
      )
      for (i in seq_len(m)) {
        preceding <- if (i == 1L) {
          0
        } else {
          colSums(
            first_magnitude[seq_len(i - 1L), , drop = FALSE] *
              abs_R[seq_len(i - 1L), i]
          )
        }
        first_magnitude[i, ] <-
          (abs(z_checked[i, ]) + preceding) / abs(R[i, i])
      }
      solution_magnitude <- base::matrix(
        0, nrow = m, ncol = length(columns_to_certify)
      )
      for (ii in seq_len(m)) {
        i <- m - ii + 1L
        following <- if (i == m) {
          0
        } else {
          colSums(
            solution_magnitude[(i + 1L):m, , drop = FALSE] *
              abs_R[i, (i + 1L):m]
          )
        }
        solution_magnitude[i, ] <-
          (first_magnitude[i, ] + following) / abs(R[i, i])
      }
      if (any(!is.finite(solution_magnitude))) return(bad)

      accuracy_tolerance <- sqrt(.Machine$double.eps) / 4
      cheap_certificate <- function(lambda_lower) {
        if (!is.finite(lambda_lower) || lambda_lower <= 0) return(FALSE)
        for (j in seq_along(columns_to_certify)) {
          residual_max <- max(residual_bound[, j])
          residual_norm <- if (residual_max == 0) {
            0
          } else {
            residual_max * sqrt(sum((residual_bound[, j] / residual_max)^2))
          }
          error_bound <- residual_norm / lambda_lower
          w_bound <-
            (error_bound +
               comparison_gamma * solution_magnitude[, j]) /
            (1 - comparison_gamma)
          solution_bound <- w_bound / scale
          if (!is.finite(error_bound) || any(!is.finite(w_bound)) ||
              any(!is.finite(solution_bound)) ||
              any(w_bound >
                    accuracy_tolerance * pmax(1, abs(w_checked[, j]))) ||
              any(solution_bound >
                    accuracy_tolerance *
                      pmax(1, abs(solution_checked[, j])))) {
            return(FALSE)
          }
        }
        TRUE
      }

      lambda_lower <- lambda_min - max(covariance_discrepancy)
      cheap_certified <- cheap_certificate(lambda_lower)
      if (!cheap_certified && condition_is_bound == 1) {
        exact <- exact_root_condition(
          R, covariance_discrepancy, condition_tolerance
        )
        if (is.null(exact)) return(bad)
        rho_K <- exact$rho
        lambda_min <- exact$lambda_min
        condition_is_bound <- 0
        relative_error_proxy <- pmax(eta, roundoff_floor) / rho_K
        if (any(!is.finite(relative_error_proxy)) ||
            any(relative_error_proxy >= 1)) return(bad)
        lambda_lower <- lambda_min - max(covariance_discrepancy)
        cheap_certified <- cheap_certificate(lambda_lower)
      }

      if (!cheap_certified) {
        # Bound the error componentwise from the residual of the completed
        # equilibrated covariance.  The gamma envelope accounts for rounding
        # in the matrix-vector product even when its observed residual is zero.
        # This catches a small lost solution component without charging a
        # worst-case triangular recurrence to every ordinary solve.
        product_length <- max(1L, m + 1L)
        product_gamma <-
          product_length * .Machine$double.eps /
          (1 - product_length * .Machine$double.eps)
        E_inverse <- tryCatch(
          base::chol2inv(R),
          error = function(e) NULL
        )
        if (is.null(E_inverse) || any(!is.finite(E_inverse)) ||
            !is.finite(product_gamma)) return(bad)

        inverse_error_proxy <-
          8 * max(1L, m) * .Machine$double.eps / rho_K
        if (!is.finite(inverse_error_proxy) || inverse_error_proxy >= 1) {
          return(bad)
        }
        completed_residual <- abs(z_checked - E %*% w_checked) +
          covariance_discrepancy %o%
            apply(abs(w_checked), 2L, max) +
          product_gamma * (
            abs(z_checked) + abs(E) %*% abs(w_checked)
        )
        if (any(!is.finite(completed_residual))) return(bad)

        # `inverse_error_proxy` is normwise, so account separately for the
        # uncertainty in the computed inverse rather than treating it as an
        # elementwise relative bound (which would not protect entries rounded
        # near zero).
        inverse_norm_error <-
          inverse_error_proxy / (1 - inverse_error_proxy) *
          max(rowSums(abs(E_inverse)))
        component_bound <- abs(E_inverse) %*% completed_residual +
          inverse_norm_error * rep(
            colSums(completed_residual), each = m
          )
        returned_bound <- component_bound / scale
        if (any(!is.finite(component_bound)) ||
            any(!is.finite(returned_bound)) ||
            any(component_bound >
                  accuracy_tolerance * pmax(1, abs(w_checked))) ||
            any(returned_bound >
                  accuracy_tolerance * pmax(1, abs(solution_checked)))) {
          return(bad)
        }
      }
    }
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

  RTMB::ADjoint(forward, reverse, name = "K2_solution_guard")
})

#' @keywords internal
.K2_solution_guard <- local({
  atomic <- .K2_solution_guard_atomic
  function(solution, R, w, z, scale, covariance_discrepancy,
           rho, lambda_min, condition_is_bound) {
    m <- nrow(R)
    nrhs <- if (is.null(dim(solution))) 1L else ncol(solution)
    solution_vec <- as.vector(solution)
    payload <- RTMB::AD(numeric(
      5L + m * m + 3L * length(solution_vec) + 2L * m
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
    cursor <- cursor + length(solution_vec)
    payload[cursor:(cursor + m - 1L)] <- as.vector(scale)
    cursor <- cursor + m
    payload[cursor:(cursor + m - 1L)] <- covariance_discrepancy
    cursor <- cursor + m
    payload[cursor:(cursor + 2L)] <- c(
      rho, lambda_min, condition_is_bound
    )
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
  covariance_is_ad <-
    inherits(packed, "advector") || inherits(packed, "adsparse") ||
    (!is.null(S) &&
      (inherits(S, "advector") || inherits(S, "adsparse")))
  factorization <- .weighted_gram_chol(packed, S)
  factorization$covariance_is_ad <- covariance_is_ad

  if (!inherits(factorization$R, "advector") &&
      any(!is.finite(c(factorization$R, factorization$scale)))) {
    stop(
      paste(
        "K2 is singular, indefinite, or numerically ill-conditioned",
        "at this evaluation."
      ),
      call. = FALSE
    )
  }
  factorization
}

#' @keywords internal
.K2_factor_solve <- function(K2_factor, tvec, parameter_vector, rhs) {
  factorization <- .K2_factor_chol(K2_factor, tvec, parameter_vector)
  terminal <- .K2_terminal_root_guard(
    factorization$R, factorization$covariance_discrepancy,
    factorization$covariance_is_ad
  )
  R <- terminal$R
  if (!inherits(R, "advector") &&
      any(!is.finite(c(
        R, terminal$rho, terminal$lambda_min,
        terminal$condition_is_bound
      )))) {
    stop("K2 solve would lose all numerical accuracy.", call. = FALSE)
  }
  z <- rhs / factorization$scale
  w <- solve(R, solve(t(R), z))
  solution <- w / factorization$scale
  solution <- .K2_solution_guard(
    solution, R, w, z, factorization$scale,
    factorization$covariance_discrepancy,
    terminal$rho, terminal$lambda_min,
    terminal$condition_is_bound
  )
  if (!inherits(solution, "advector") && any(!is.finite(solution))) {
    stop("K2 solve failed numerical accuracy certification.", call. = FALSE)
  }
  solution
}

#' @keywords internal
.K2_factor_logdet <- function(K2_factor, tvec, parameter_vector) {
  factorization <- .K2_factor_chol(K2_factor, tvec, parameter_vector)
  terminal <- .K2_terminal_root_guard(
    factorization$R, factorization$covariance_discrepancy,
    factorization$covariance_is_ad
  )
  R <- terminal$R
  if (!inherits(R, "advector") && any(!is.finite(R))) {
    stop(
      "K2 log-determinant derivatives would lose all numerical accuracy.",
      call. = FALSE
    )
  }
  2 * sum(log(diag(R))) +
    2 * sum(log(factorization$scale))
}
