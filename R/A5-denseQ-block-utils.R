# Internal helpers for exact dense-Q and factored K3K3 operators on wrappers
# built from independent blocks.
# used in R/A-CGF_IIDReplicatesCGF.R and R/D-ConcatenationCGF.R



# Denote the 3rd- and 4th-derivative tensors of the CGF as
# K_{ijk}(t)   = partial^3 K(t) / partialt_i partialt_j partialt_k
# K_{ijkl}(t)  = partial^4 K(t) / partialt_i partialt_j ∂t_k partialt_l

# Then for a matrix Q, the 3 higher-order contractions are:
#   K4_AABB(Q)
#   = \sum_{a,b,c,d} K_{abcd} Q_{ab} Q_{cd}
#   K3K3_AABBCC(Q)
#   = sum_{a,b,c,i,j,k} K_{abc} K_{ijk} Q_{ab} Q_{ij} Q_{ck}
#   K3K3_ABCABC(Q)
#   = sum_{a,b,c,i,j,k} K_{abc} K_{ijk} Q_{ai} Q_{bj} Q_{ck}

# Now suppose we split a random vector into independent blocks: X = (X^(1), ..., X^(B)),
# then the overall 3rd-cumulant tensor is block-sparse:
# *** K_{abc} = 0 unless a,b,c all belong to the same block
# and similarly
# *** K_{abcd} = 0 unless a,b,c,d all belong to the same block

      ### therefore K4_AABB is easy because from this:
      ### K4_AABB(Q) = sum_{a,b,c,d} K_{abcd} Q_{ab} Q_{cd}, we end up with
      ### K4_AABB(Q) = sum_b K4_AABB^(b)(Q_bb) ==> using only Q_bb

# K3K3_AABBCC is not straighforward
# starting from:
# K3K3_AABBCC(Q) = sum_{a,b,c,i,j,k} K_{abc} K_{ijk} Q_{ab} Q_{ij} Q_{ck}
# and with block independence,
# - a,b,c must be in one block, say block b
# - i,j,k must be in one block, say block β

# So we proceed as follows:
# sum_{b,β}
# sum_{a,a',c in block b}
#   sum_{i,j,k in block β}
#   T_b[a,a',c] T_β[i,j,k]
# (Q_bb)[a,a'] (Q_ββ)[i,j] Q_{bβ}[c,k]
     #####.... T here has nothing to do with func_T, it just used because of the word tensor
# then for each block b, we define a vector
# u_c^(b) = \sum_{a,a'} T_b[a,a',c] (Q_bb)[a,a']
# essentially: contract the first 2 indices of the block K3 tensor with the block-diagonal part Q_bb

# we end up with this:
# K3K3_AABBCC(Q)
# = \sum_{b,β} \sum_{c,k} u_c^(b) Q_{bβ}[c,k] u_k^(β)
# = \sum_{b,β} (u^(b))^T Q_{bβ} u^(β)

####### END K3K3_AABBCC


# K3K3_ABCABC is different
# starting from
# K3K3_ABCABC(Q) = \sum_{a,b,c,i,j,k} K_{abc} K_{ijk} Q_{ai} Q_{bj} Q_{ck},
# we write
# K3K3_ABCABC(Q)
# = \sum_{b,β} \sum T_b[a,b',c]
# T_β[i,j,k] Q_{bβ}[a,i] Q_{bβ}[b',j] Q_{bβ}[c,k]

# using  M = Q_{bβ}, one block-pair contribution has the form
# pair(b,β) = \sum T_b[a,b',c] T_β[i,j,k] M[a,i] M[b',j] M[c,k]
# and K3K3_ABCABC(Q) = sum_{b,β} pair(b,β)
#### For dense Q, ABCABC cannot be reduced to one vector per block in the same
#### clean way as AABBCC.

####### END K3K3_ABCABC

# For the factored forms, if Q = A diag(d) A^T,
# then the block (b,β) of Q is
# Q_{bβ} = A_b diag(d) A_β^T, where A_b ==> the rows of A belonging to block b
# ==> we can use the dense formulas; we just compute each needed block from the factorization.
# For AABBCC:
# 1) build each u^(b) using Q_bb = A_b diag(d) A_b^T
# 2) stack them into one global u
# 3) AABBCC = u^T Q u
#           = u^T A diag(d) A^T u
# we can define z = A^T u, then AABBCC = z^T diag(d) z = \sum_r d_r z_r^2
# For ABCABC, use M_{bβ} = A_b diag(d) A_β^T inside the block-pair
# formula.



##### Other related things:
# The base formulas e.g AABBCC as Σ_{m2} d[m2] (sum_{m1} d[m1] K3(A[,m1], A[,m1], A[,m2]))^2,
# if we expand this and consider that for an independent-blocks, K3(...) only gets contributions
# when all three vector arguments hit the same block, we'd get the formula
# AABBCC = sum_{b,β} (u^(b))^T Q_{bβ} u^(β)  ....### I will need to be sure of this???

# So are the base methods better? I do not think so!
# The wrappers we define here work with block size d, not full dimension length(tvec)=N ... see in R/A-CGF_IIDReplicatesCGF.R
# SO this should be better for
# *** many independent blocks
# *** small/moderate block dimension
# *** dense Q

# But again .extract_symmetric_K3_tensor() here has 3 nested loops over block dimension d.
# Even with symmetry, it is still basically O(d^3) K3operator evaluations per block.
# So it is definitely not computationally the perfect; but i think it may still be
# acceptable here because it is over the block size, not the full global dimension.

# I'd expect that in iidReplicates or concatenation:
# *** global dimension may be large
# *** but each independent block should be small (not necessarily)
# SO we'd be extracting K3 once per block: B * O(d^3), then block contractions
# Should still be better than generic base route, which ignores independence

##### Most obvious improvements:
# We do not reconstruct the full block K3 tensor if the child already provides enough structure
# For example:
# - if the child has its own fast K3K3operatorAABBCC_factored/K3K3operatorABCABC_factored, we'd use that where possible
# - if the block CGF is diagonal/vectorized/..., contract directly without building all d^3 entries
# - another path, cache extracted block tensors when the same block CGF and tvec[idx] are reused within one call path???? maybe not??




.ad_type_scale <- function(param) {
  if (length(param) > 0) param[1] else 0
}

.ad_zero_scalar <- function(param) {
  0 * .ad_type_scale(param)
}

.ad_zero_vector <- function(n, param) {
  numeric(n) * .ad_type_scale(param)
}

.ad_zero_array <- function(dim, param) {
  array(0, dim = dim) * .ad_type_scale(param)
}

# Re-express Q = A diag(dvec) A' with columns on a comparable scale.  The
# regularization keeps the transformation differentiable when a column of A is
# exactly zero; it is not a covariance jitter because the transformed factors
# represent exactly the same Q:
#
#   C_j = A_j / s_j,  w_j = d_j s_j^2,
#   C diag(w) C' = A diag(d) A'.
#
# Keeping dvec separate (rather than multiplying A by sqrt(dvec)) is important:
# values, gradients and Hessians then remain polynomial at zero or underflowed
# weights.  The scale is evaluated on the tape, so a parameter-dependent factor
# remains balanced when the tape is evaluated away from its recording point.
.balance_factored_Q <- function(A, dvec, where = "factored contraction") {
  A_dim <- dim(A)
  if (length(A_dim) != 2L) {
    stop(where, ": A must be matrix-like.", call. = FALSE)
  }
  if (A_dim[2L] != length(dvec)) {
    stop(
      where,
      ": Column/weight mismatch: ncol(A) must equal length(dvec).",
      call. = FALSE
    )
  }

  r <- length(dvec)
  if (r == 0L) return(list(A = A, d = dvec))

  # Use a mean absolute column scale rather than an RMS scale: squaring the raw
  # factor can overflow even when the represented covariance is ordinary
  # because dvec supplies the reciprocal scale.  A half-mean leaves a factor-
  # two overflow margin while keeping every finite column on a useful scale.
  # The small positive floor gives finite arithmetic at an exactly zero
  # column; it is not added to Q.
  n <- A_dim[1L]
  scales <- as.vector(colSums(abs(A) / (2 * max(1L, n)))) +
    sqrt(.Machine$double.xmin)
  balanced_A <- t(t(A) / scales)
  balanced_d <- (dvec * scales) * scales

  list(A = balanced_A, d = balanced_d)
}


.symmetric_K3_count <- function(d) {
  as.double(d) * (d + 1) * (d + 2) / 6
}

.use_direct_factored_rank <- function(block_dims, r,
                                      method = c("AABBCC", "ABCABC")) {
  method <- match.arg(method)
  if (r == 0L || length(block_dims) == 0L) return(FALSE)

  # One direction costs exactly one child call per block for either route and
  # avoids rebuilding a cancellation-prone coordinate projection.
  if (length(block_dims) > 1L && r == 1L) return(TRUE)

  # This route is intended for genuinely thin factors.  Near-square factors
  # have smaller derivative tapes through coordinate extraction even when the
  # raw callback counts are similar.
  if (r > floor(min(block_dims) / 2)) return(FALSE)

  coordinate_callbacks <- sum(.symmetric_K3_count(block_dims))
  rank_callbacks <- length(block_dims) * if (method == "AABBCC") {
    as.double(r) * r
  } else {
    .symmetric_K3_count(r)
  }

  rank_callbacks < coordinate_callbacks
}


# Exact rank-space contractions for an independent-block third cumulant.  The
# block contributions are accumulated before squaring, as required by the
# tensor identity and for cancellation stability.
.block_K3K3_AABBCC_rank <- function(K3fun_list, tvec_blocks, row_blocks,
                                    dvec, param) {
  B <- length(row_blocks)
  r <- length(dvec)
  if (length(K3fun_list) != B || length(tvec_blocks) != B) {
    stop("Internal block-list mismatch in AABBCC rank contraction.")
  }

  total <- .ad_zero_scalar(param)
  for (m2 in seq_len(r)) {
    inner <- .ad_zero_scalar(param)
    for (m1 in seq_len(r)) {
      across_blocks <- .ad_zero_scalar(param)
      for (b in seq_len(B)) {
        A_b <- row_blocks[[b]]
        across_blocks <- across_blocks + K3fun_list[[b]](
          tvec_blocks[[b]], param,
          as.vector(A_b[, m1]),
          as.vector(A_b[, m1]),
          as.vector(A_b[, m2])
        )
      }
      inner <- inner + dvec[m1] * across_blocks
    }
    total <- total + (dvec[m2] * inner) * inner
  }

  total
}

.block_K3K3_ABCABC_rank <- function(K3fun_list, tvec_blocks, row_blocks,
                                    dvec, param) {
  B <- length(row_blocks)
  r <- length(dvec)
  if (length(K3fun_list) != B || length(tvec_blocks) != B) {
    stop("Internal block-list mismatch in ABCABC rank contraction.")
  }

  total <- .ad_zero_scalar(param)
  for (m1 in seq_len(r)) {
    for (m2 in m1:r) {
      for (m3 in m2:r) {
        across_blocks <- .ad_zero_scalar(param)
        for (b in seq_len(B)) {
          A_b <- row_blocks[[b]]
          across_blocks <- across_blocks + K3fun_list[[b]](
            tvec_blocks[[b]], param,
            as.vector(A_b[, m1]),
            as.vector(A_b[, m2]),
            as.vector(A_b[, m3])
          )
        }
        multiplicity <- if (m1 == m3) {
          1
        } else if (m1 == m2 || m2 == m3) {
          3
        } else {
          6
        }
        weighted_square <-
          ((dvec[m1] * across_blocks) * dvec[m2]) *
            (dvec[m3] * across_blocks)
        total <- total + multiplicity * weighted_square
      }
    }
  }

  total
}

.extract_symmetric_K3_tensor <- function(K3fun,
                                         tvec,
                                         param,
                                         block_dim,
                                         basis = NULL) {
  if (block_dim < 1L) stop("'block_dim' must be >= 1.")
  if (is.null(basis)) basis <- diag(1, block_dim)

  T3 <- .ad_zero_array(c(block_dim, block_dim, block_dim), param)

  for (i in seq_len(block_dim)) {
    ei <- basis[, i]
    for (j in i:block_dim) {
      ej <- basis[, j]
      for (k in j:block_dim) {
        ek <- basis[, k]
        val <- K3fun(tvec, param, ei, ej, ek)

        perms <- unique(rbind(
          c(i, j, k),
          c(i, k, j),
          c(j, i, k),
          c(j, k, i),
          c(k, i, j),
          c(k, j, i)
        ))

        for (r in seq_len(nrow(perms))) {
          T3[perms[r, 1L], perms[r, 2L], perms[r, 3L]] <- val
        }
      }
    }
  }

  T3
}

# this builds T[:,:,k]
.extract_K3_slices <- function(K3fun,
                               tvec,
                               param,
                               block_dim,
                               basis = NULL) {
  T3 <- .extract_symmetric_K3_tensor(
    K3fun = K3fun,
    tvec = tvec,
    param = param,
    block_dim = block_dim,
    basis = basis
  )

  lapply(seq_len(block_dim), function(k) {
    matrix(T3[, , k], nrow = block_dim, ncol = block_dim)
  })
}


# computes u_k = sum_{a,a'} T[a,a',k] Q_block[a,a']
.k3_slices_to_aabbcc_vector <- function(k3_slices, Q_block, param) {
  d <- length(k3_slices)
  if (nrow(Q_block) != d || ncol(Q_block) != d) {
    stop("Q_block dimension mismatch in .k3_slices_to_aabbcc_vector().")
  }

  out <- .ad_zero_vector(d, param)
  for (k in seq_len(d)) {
    out[k] <- sum(k3_slices[[k]] * Q_block)
  }
  out
}

# computes one \sum T_left T_right M M M   block-pair term
.k3_slices_abcabc_pair <- function(k3_slices_left,
                                   k3_slices_right,
                                   M,
                                   param) {
  d_left <- length(k3_slices_left)
  d_right <- length(k3_slices_right)

  if (nrow(M) != d_left || ncol(M) != d_right) {
    stop("Block dimension mismatch in .k3_slices_abcabc_pair().")
  }

  total <- .ad_zero_scalar(param)
  Mt <- t(M)

  for (k in seq_len(d_right)) {
    S <- M %*% k3_slices_right[[k]] %*% Mt
    for (c in seq_len(d_left)) {
      total <- total + M[c, k] * sum(k3_slices_left[[c]] * S)
    }
  }

  total
}

.k3_slices_abcabc_from_dense_Q <- function(k3_by_block,
                                           block_indices,
                                           Q,
                                           param) {
  B <- length(block_indices)
  total <- .ad_zero_scalar(param)

  for (i in seq_len(B)) {
    idx_i <- block_indices[[i]]
    for (j in i:B) {
      idx_j <- block_indices[[j]]
      term <- .k3_slices_abcabc_pair(
        k3_slices_left = k3_by_block[[i]],
        k3_slices_right = k3_by_block[[j]],
        M = Q[idx_i, idx_j, drop = FALSE],
        param = param
      )
      total <- total + if (i == j) term else 2 * term
    }
  }

  total
}

.factor_block_matrix <- function(B_left, dvec, B_right) {
  if (ncol(B_left) != length(dvec) || ncol(B_right) != length(dvec)) {
    stop("Column/weight mismatch in .factor_block_matrix().")
  }

  scaled_left <- t(t(B_left) * dvec)
  scaled_left %*% t(B_right)
}

.k3_slices_abcabc_from_factored_Q <- function(k3_by_block,
                                              row_blocks,
                                              dvec,
                                              param) {
  B <- length(row_blocks)
  r <- length(dvec)

  factor_columns <- vapply(row_blocks, ncol, integer(1L))
  if (any(factor_columns != r)) {
    stop("Column/weight mismatch in .factor_block_matrix().")
  }
  if (B == 0L || r == 0L) return(.ad_zero_scalar(param))

  block_dims <- as.integer(vapply(row_blocks, nrow, integer(1L)))
  balanced <- .balance_factored_Q(
    do.call(rbind, row_blocks),
    dvec,
    ".k3_slices_abcabc_from_factored_Q"
  )
  ends <- cumsum(block_dims)
  starts <- ends - block_dims + 1L
  row_blocks <- Map(
    function(first, last) balanced$A[first:last, , drop = FALSE],
    starts,
    ends
  )
  dvec <- balanced$d

  # Once the block tensors have been extracted, either expand all block pairs
  # or aggregate the tensor in factor rank before the final square.  The latter
  # is linear in the number of blocks and is essential for a thin factor over
  # many scalar/small blocks.
  block_dims_double <- as.double(block_dims)
  n_factor_triples <- .symmetric_K3_count(r)
  rank_cost <- n_factor_triples + sum(
    block_dims_double^3 * r +
      block_dims_double^2 * r^2 +
      block_dims_double * n_factor_triples
  )

  suffix_d1 <- rev(cumsum(rev(block_dims_double)))
  suffix_d2 <- rev(cumsum(rev(block_dims_double^2)))
  suffix_d3 <- rev(cumsum(rev(block_dims_double^3)))
  n_remaining <- rev(seq_len(B))
  pair_cost <- sum(
    block_dims_double * r * n_remaining +
      block_dims_double * r * suffix_d1 +
      block_dims_double * suffix_d3 +
      block_dims_double^2 * suffix_d2 +
      block_dims_double^3 * suffix_d1
  )

  common_workspace <-
    sum(block_dims_double^3) + sum(block_dims_double) * r
  rank_workspace <- common_workspace +
    4 * n_factor_triples + as.double(r)^2 + max(block_dims_double) * r
  pair_workspace <- common_workspace +
    2 * max(block_dims_double)^2 + max(block_dims_double) * r

  # Reject a cheaper arithmetic route if its explicit index/accumulator
  # workspace would exceed three times the block-pair alternative.
  # Very small ranks are also kept in rank space even when they are full rank.
  # In that bounded case the block-pair formula can round several large pair
  # contractions independently before they cancel, while rank aggregation
  # performs the cancellation once.  The absolute r <= 3 guard cannot affect
  # the ordinary moderate/large full-rank path.
  bounded_full_rank <- B > 1L && r <= 3L
  use_rank <- n_factor_triples <= .Machine$integer.max && (
    bounded_full_rank || (
      r <= floor(sum(block_dims) / 2) &&
        rank_cost < pair_cost &&
        rank_workspace <= 3 * max(1, pair_workspace)
    )
  )

  if (use_rank) {
    n_factor_triples <- as.integer(n_factor_triples)
    triple_m1 <- integer(n_factor_triples)
    triple_m2 <- integer(n_factor_triples)
    triple_m3 <- integer(n_factor_triples)
    multiplicity <- integer(n_factor_triples)
    q <- 0L
    for (m1 in seq_len(r)) {
      for (m2 in m1:r) {
        for (m3 in m2:r) {
          q <- q + 1L
          triple_m1[q] <- m1
          triple_m2[q] <- m2
          triple_m3[q] <- m3
          multiplicity[q] <- if (m1 == m3) {
            1L
          } else if (m1 == m2 || m2 == m3) {
            3L
          } else {
            6L
          }
        }
      }
    }

    aggregate_values <- lapply(
      seq_len(n_factor_triples),
      function(q) .ad_zero_scalar(param)
    )
    for (i in seq_len(B)) {
      A_i <- row_blocks[[i]]
      for (k in seq_len(nrow(A_i))) {
        projected_slice <- crossprod(
          A_i,
          k3_by_block[[i]][[k]] %*% A_i
        )
        for (q in seq_len(n_factor_triples)) {
          value <-
            projected_slice[triple_m1[q], triple_m2[q]] *
            A_i[k, triple_m3[q]]
          aggregate_values[[q]] <- aggregate_values[[q]] + value
        }
      }
    }

    total <- .ad_zero_scalar(param)
    for (q in seq_len(n_factor_triples)) {
      m1 <- triple_m1[q]
      m2 <- triple_m2[q]
      m3 <- triple_m3[q]
      value <- aggregate_values[[q]]
      weighted_square <-
        ((dvec[m1] * value) * dvec[m2]) * (dvec[m3] * value)
      total <- total + multiplicity[q] * weighted_square
    }
    return(total)
  }

  total <- .ad_zero_scalar(param)

  for (i in seq_len(B)) {
    A_i <- row_blocks[[i]]
    for (j in i:B) {
      M <- .factor_block_matrix(A_i, dvec, row_blocks[[j]])
      term <- .k3_slices_abcabc_pair(
        k3_slices_left = k3_by_block[[i]],
        k3_slices_right = k3_by_block[[j]],
        M = M,
        param = param
      )
      total <- total + if (i == j) term else 2 * term
    }
  }

  total
}
