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
#### ABCABC cannot be reduced to one vector per block in the same clean way as AABBCC.

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
# For ABCABC, we just use M_{bβ} = A_b diag(d) A_β^T inside the same block-pair formula



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
