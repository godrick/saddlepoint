# R/ConcatenationCGF.R
#
# A CGF wrapper for concatenating *independent* random vectors whose CGFs act on
# disjoint sub-blocks of tvec:
#
#   tvec = (t^(1), ..., t^(L))  with dim(t^(i)) = component_dims[i]
#   K(tvec;theta) = sum_i K_i(t^(i);theta)
#
# This differs from sumOfIndependentCGF(): in sumOfIndependentCGF all children
# share the SAME tvec (same dimension) and are added. Here each child sees its
# OWN slice of tvec (block concatenation).
#
# The returned CGF is wrapped with iidReplicatesCGF using
#   block_size = sum(component_dims).
# With iidReps = "any" (default), the number of i.i.d. blocks is inferred from
# length(tvec) (must be a multiple of block_size).







.concatenationCGF_internal <- function(cgf_list,
                                       component_dims,
                                       ...){

  # checks
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if (any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE))) stop("Every element of 'cgf_list' must be an object of class 'CGF'.")


  L <- length(cgf_list)
  dims <- as.integer(component_dims)
  if (length(dims) != L) stop("Internal error: component_dims length mismatch.")

  total_dim <- sum(dims)

  # Precompute indices once
  ends   <- cumsum(dims)
  starts <- ends - dims + 1L
  idx_list <- Map(seq.int, starts, ends)







  K_list      <- lapply(cgf_list, function(cg) cg$K)
  K1_list     <- lapply(cgf_list, function(cg) cg$K1)
  K2_list     <- lapply(cgf_list, function(cg) cg$K2)
  K3op_list   <- lapply(cgf_list, function(cg) cg$K3operator)
  K4op_list   <- lapply(cgf_list, function(cg) cg$K4operator)

  K2op_list       <- lapply(cgf_list, function(cg) cg$K2operator)
  K2opAK2AT_list  <- lapply(cgf_list, function(cg) cg$K2operatorAK2AT)
  K2_factor_list  <- lapply(cgf_list, .K2_factor_method)

  K2_solve_list   <- lapply(cgf_list, function(cg) cg$K2_solve)
  logdet_list     <- lapply(cgf_list, function(cg) cg$logdetK2)

  K4AABB_list     <- lapply(cgf_list, function(cg) cg$K4operatorAABB)
  K3K3AABBCC_list <- lapply(cgf_list, function(cg) cg$K3K3operatorAABBCC)
  K3K3ABCABC_list <- lapply(cgf_list, function(cg) cg$K3K3operatorABCABC)
  K4AABB_factored_list <- lapply(cgf_list, function(cg) cg$.private_api$K4operatorAABB_factored)
  K3K3AABBCC_factored_list <- lapply(
    cgf_list,
    function(cg) cg$.private_api$K3K3operatorAABBCC_factored
  )
  K3K3ABCABC_factored_list <- lapply(
    cgf_list,
    function(cg) cg$.private_api$K3K3operatorABCABC_factored
  )
  K4AABB_delegate_safe <- vapply(
    K4AABB_factored_list,
    .factored_delegate_is_safe,
    logical(1)
  )
  K3K3AABBCC_delegate_safe <- vapply(
    K3K3AABBCC_factored_list,
    .factored_delegate_is_safe,
    logical(1)
  )
  K3K3ABCABC_delegate_safe <- vapply(
    K3K3ABCABC_factored_list,
    .factored_delegate_is_safe,
    logical(1)
  )

  ineq_list       <- lapply(cgf_list, function(cg) cg$ineq_constraint)

  #
  tilting_list <- lapply(cgf_list, function(cg) cg$.private_api$tilting_exponent)
  negll_list   <- lapply(cgf_list, function(cg) cg$.private_api$neg_ll)
  funcT_list   <- lapply(cgf_list, function(cg) cg$.private_api$func_T)

  # Analytic t-hat: only if ALL children have it
  has_analytic_vec <- vapply(cgf_list, function(cg) isTRUE(cg$has_analytic_tvec_hat), logical(1))
  analytic_hat_list <- if (all(has_analytic_vec)) {
    lapply(cgf_list, function(cg) cg$.private_api$analytic_tvec_hat_func)
  } else {
    NULL
  }
  basis_list <- lapply(dims, function(d) diag(1, d))

  validate_factored_K3K3 <- function(tvec, A, dvec, where) {
    if (length(tvec) != total_dim) {
      stop(where, ": tvec length mismatch.")
    }
    A_dim <- dim(A)
    if (length(A_dim) != 2L || A_dim[1L] != total_dim) {
      stop(
        where,
        ": A must be matrix-like with nrow(A) == length(tvec)."
      )
    }
    if (A_dim[2L] != length(dvec)) {
      stop(
        where,
        ": Column/weight mismatch: ncol(A) must equal length(dvec)."
      )
    }
    length(dvec)
  }












  # # K => sum of child K, but each sees sub-block
  # Kfun <- function(tvec, param) {
  #   # param is global, each CGF in cgf_list knows what to do with it
  #   # tvec must have length = total_dim
  #   if (length(tvec)!= total_dim) stop(sprintf("`tvec` has length %d, expected %d from component_dims sum.", length(tvec), total_dim))
  #   total <- 0
  #   current_start <- 1
  #   for (i in seq_along(cgf_list)) {
  #     len_i <- component_dims[i]
  #     t_sub <- tvec[current_start:(current_start+len_i-1)]
  #     total <- total + cgf_list[[i]]$K(t_sub, param)
  #     current_start <- current_start + len_i
  #   }
  #   total
  # }

  K <- function(tvec, param) {
    if (length(tvec) != total_dim) {
      stop(sprintf("`tvec` has length %d, expected %d (= sum(component_dims)).",
                   length(tvec), total_dim))
    }
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K_list[[i]](tvec[idx], param)
    }
    total
  }


  K1 <- function(tvec, param) {
    if (length(tvec) != total_dim) {
      stop(sprintf("`tvec` has length %d, expected %d (= sum(component_dims)).",
                   length(tvec), total_dim))
    }
    out <- numeric(total_dim) * param[1]
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      out[idx] <- K1_list[[i]](tvec[idx], param)
    }
    out
  }

  K2 <- function(tvec, param) {
    if (length(tvec) != total_dim) {
      stop(sprintf("`tvec` has length %d, expected %d (= sum(component_dims)).",
                   length(tvec), total_dim))
    }
    accum <- matrix(0, nrow = total_dim, ncol = total_dim) * param[1]
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      d_i <- dims[i]
      k2_sub <- K2_list[[i]](tvec[idx], param)

      # univariate return types (scalar -> 1x1)
      if (is.null(dim(k2_sub))) k2_sub <- matrix(k2_sub, nrow = d_i, ncol = d_i)
      accum[idx, idx] <- k2_sub
    }
    accum
  }






  # # K3operator => sum of child K3operator
  # K3opfun <- function(tvec, param, v1, v2, v3) {
  #   if (length(tvec)!= total_dim || length(v1)!= total_dim ||
  #       length(v2)!= total_dim || length(v3)!= total_dim) {
  #     stop("dimension mismatch in K3operator arguments.")
  #   }
  #   total <- 0
  #   current_start <- 1
  #   for (i in seq_along(cgf_list)) {
  #     len_i <- component_dims[i]
  #     idx <- current_start:(current_start+len_i-1)
  #     total <- total + cgf_list[[i]]$K3operator(tvec[idx], param, v1[idx], v2[idx], v3[idx])
  #     current_start <- current_start + len_i
  #   }
  #   total
  # }
  K3operator <- function(tvec, param, v1, v2, v3) {
    if (length(tvec) != total_dim ||
        length(v1)   != total_dim ||
        length(v2)   != total_dim ||
        length(v3)   != total_dim) {
      stop("dimension mismatch in K3operator arguments.")
    }
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K3op_list[[i]](tvec[idx], param, v1[idx], v2[idx], v3[idx])
    }
    total
  }









  # # K4operator => sum across sub-block operators
  # K4opfun <- function(tvec, param, v1, v2, v3, v4) {
  #   if (length(tvec)!= total_dim || length(v1)!= total_dim ||
  #       length(v2)!= total_dim || length(v3)!= total_dim || length(v4)!= total_dim) {
  #     stop("dimension mismatch in K4operator arguments.")
  #   }
  #   total <- 0
  #   current_start <- 1
  #   for (i in seq_along(cgf_list)) {
  #     len_i <- component_dims[i]
  #     idx <- current_start:(current_start+len_i-1)
  #     total <- total + cgf_list[[i]]$K4operator(tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
  #     current_start <- current_start + len_i
  #   }
  #   total
  # }

  K4operator <- function(tvec, param, v1, v2, v3, v4) {
    if (length(tvec) != total_dim ||
        length(v1)   != total_dim ||
        length(v2)   != total_dim ||
        length(v3)   != total_dim ||
        length(v4)   != total_dim) {
      stop("dimension mismatch in K4operator arguments.")
    }
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K4op_list[[i]](tvec[idx], param, v1[idx], v2[idx], v3[idx], v4[idx])
    }
    total
  }



  # ---------------------------------------------------------------------------
  # Derived methods (block structure)
  # ---------------------------------------------------------------------------
  tilting_exponent <- function(tvec, param) {
    if (length(tvec) != total_dim) stop("tilting_exponent: tvec length mismatch.")
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + tilting_list[[i]](tvec[idx], param)
    }
    total
  }

  # func_T: sum of child func_T (blockwise).
  func_T <- function(tvec, param) {
    if (length(tvec) != total_dim) stop("func_T: tvec length mismatch.")
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + funcT_list[[i]](tvec[idx], param)
    }
    total
  }

  neg_ll <- function(tvec, param) {
    if (length(tvec) != total_dim) stop("neg_ll: tvec length mismatch.")
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + negll_list[[i]](tvec[idx], param)
    }
    total
  }


  K2operator <- function(tvec, param, x, y) {
    if (length(tvec) != total_dim || length(x) != total_dim || length(y) != total_dim) {
      stop("K2operator: dimension mismatch.")
    }
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K2op_list[[i]](tvec[idx], param, x[idx], y[idx])
    }
    total
  }


  K2operatorAK2AT <- function(tvec, param, Bmat) {
    if (length(tvec) != total_dim) stop("K2operatorAK2AT: tvec length mismatch.")
    if (ncol(Bmat) != total_dim) {
      stop("K2operatorAK2AT: Bmat must have ncol == length(tvec). ",
           "Got ncol(Bmat)=", ncol(Bmat), ", length(tvec)=", total_dim, ".")
    }
    r <- nrow(Bmat)
    out <- matrix(0, nrow = r, ncol = r) * param[1]
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      out <- out + K2opAK2AT_list[[i]](tvec[idx], param, Bmat[, idx, drop = FALSE])
    }
    out
  }

  K2_factor <- NULL
  if (any(vapply(K2_factor_list, is.function, logical(1)))) {
    K2_factor <- function(tvec, param, Bmat) {
      if (length(tvec) != total_dim || ncol(Bmat) != total_dim) {
        stop("K2_factor: dimension mismatch.")
      }
      terms <- list()
      for (i in seq_len(L)) {
        idx <- idx_list[[i]]
        child_terms <- if (is.function(K2_factor_list[[i]])) {
          K2_factor_list[[i]](tvec[idx], param, Bmat[, idx, drop = FALSE])
        } else {
          .K2_dense_term(
            K2opAK2AT_list[[i]](tvec[idx], param, Bmat[, idx, drop = FALSE])
          )
        }
        terms[[length(terms) + 1L]] <- child_terms
      }
      unlist(terms, recursive = FALSE)
    }
  }


  K2_solve <- function(tvec, param, rhs) {
    if (length(tvec) != total_dim) stop("K2_solve: tvec length mismatch.")

    # vector rhs
    if (is.null(dim(rhs))) {
      if (length(rhs) != total_dim) stop("K2_solve: rhs length mismatch.")
      out <- numeric(total_dim)*param[1]
      for (i in seq_len(L)) {
        idx <- idx_list[[i]]
        out[idx] <- K2_solve_list[[i]](tvec[idx], param, rhs[idx])
      }
      return(out)
    }

    # matrix rhs
    if (nrow(rhs) != total_dim) stop("K2_solve: rhs nrow mismatch.")
    k <- ncol(rhs)
    out <- matrix(0, nrow = total_dim, ncol = k)*param[1]
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      out[idx, ] <- K2_solve_list[[i]](tvec[idx], param, rhs[idx, , drop = FALSE])
    }
    out
  }

  logdetK2 <- function(tvec, param) {
    if (length(tvec) != total_dim) stop("logdetK2: tvec length mismatch.")
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + logdet_list[[i]](tvec[idx], param)
    }
    total
  }


  # ---------------------------------------------------------------------------
  # Higher-order operators for independent-block concatenation.
  # K4operatorAABB only sees diagonal sub-blocks, but the K3K3 operators also
  # require off-diagonal block couplings from Q.
  # ---------------------------------------------------------------------------

  K4operatorAABB <- function(tvec, param, Q) {
    if (length(tvec) != total_dim) stop("K4operatorAABB: tvec length mismatch.")
    total <- 0
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K4AABB_list[[i]](tvec[idx], param,
                                        Q[idx, idx, drop = FALSE])
    }
    total
  }

  K3K3operatorAABBCC <- function(tvec, param, Q) {
    if (length(tvec) != total_dim) stop("K3K3operatorAABBCC: tvec length mismatch.")
    if (L == 1L) {
      return(K3K3AABBCC_list[[1L]](tvec, param, Q))
    }
    u <- .ad_zero_vector(total_dim, param)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      k3_slices <- .extract_K3_slices(
        K3fun = K3op_list[[i]],
        tvec = tvec[idx],
        param = param,
        block_dim = dims[i],
        basis = basis_list[[i]]
      )
      u[idx] <- .k3_slices_to_aabbcc_vector(k3_slices, Q[idx, idx, drop = FALSE], param)
    }

    sum(u * as.vector(Q %*% u))
  }

  K3K3operatorABCABC <- function(tvec, param, Q) {
    if (length(tvec) != total_dim) stop("K3K3operatorABCABC: tvec length mismatch.")
    if (L == 1L) {
      return(K3K3ABCABC_list[[1L]](tvec, param, Q))
    }
    k3_by_block <- vector("list", L)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      k3_by_block[[i]] <- .extract_K3_slices(
        K3fun = K3op_list[[i]],
        tvec = tvec[idx],
        param = param,
        block_dim = dims[i],
        basis = basis_list[[i]]
      )
    }

    .k3_slices_abcabc_from_dense_Q(
      k3_by_block = k3_by_block,
      block_indices = idx_list,
      Q = Q,
      param = param
    )
  }

  K4operatorAABB_factored <- function(tvec, param, A, dvec) {
    r <- validate_factored_K3K3(
      tvec, A, dvec, "K4operatorAABB_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))
    total <- .ad_zero_scalar(param)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      total <- total + K4AABB_factored_list[[i]](
        tvec[idx], param, A[idx, , drop = FALSE], dvec
      )
    }

    total
  }
  K4operatorAABB_factored <- .factored_delegate_mark(
    K4operatorAABB_factored, all(K4AABB_delegate_safe)
  )

  K3K3operatorAABBCC_factored <- function(tvec, param, A, dvec) {
    r <- validate_factored_K3K3(
      tvec, A, dvec, "K3K3operatorAABBCC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))
    if (L == 1L && K3K3AABBCC_delegate_safe[[1L]]) {
      return(K3K3AABBCC_factored_list[[1L]](
        tvec, param, A, dvec
      ))
    }

    if (.use_direct_factored_rank(dims, r, "AABBCC")) {
      balanced <- .balance_factored_Q(
        A, dvec, "concatenation K3K3operatorAABBCC_factored"
      )
      return(.block_K3K3_AABBCC_rank(
        K3fun_list = K3op_list,
        tvec_blocks = lapply(idx_list, function(idx) tvec[idx]),
        row_blocks = lapply(
          idx_list,
          function(idx) balanced$A[idx, , drop = FALSE]
        ),
        dvec = balanced$d,
        param = param
      ))
    }

    balanced <- .balance_factored_Q(
      A, dvec, "concatenation K3K3operatorAABBCC_factored fallback"
    )
    A <- balanced$A
    dvec <- balanced$d

    u <- .ad_zero_vector(total_dim, param)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      A_block <- A[idx, , drop = FALSE]
      Qii <- .factor_block_matrix(A_block, dvec, A_block)
      k3_slices <- .extract_K3_slices(
        K3fun = K3op_list[[i]],
        tvec = tvec[idx],
        param = param,
        block_dim = dims[i],
        basis = basis_list[[i]]
      )
      u[idx] <- .k3_slices_to_aabbcc_vector(k3_slices, Qii, param)
    }

    z <- as.vector(crossprod(A, u))
    sum(dvec * z * z)
  }

  K3K3operatorABCABC_factored <- function(tvec, param, A, dvec) {
    r <- validate_factored_K3K3(
      tvec, A, dvec, "K3K3operatorABCABC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))
    if (L == 1L && K3K3ABCABC_delegate_safe[[1L]]) {
      return(K3K3ABCABC_factored_list[[1L]](
        tvec, param, A, dvec
      ))
    }

    if (.use_direct_factored_rank(dims, r, "ABCABC")) {
      balanced <- .balance_factored_Q(
        A, dvec, "concatenation K3K3operatorABCABC_factored"
      )
      return(.block_K3K3_ABCABC_rank(
        K3fun_list = K3op_list,
        tvec_blocks = lapply(idx_list, function(idx) tvec[idx]),
        row_blocks = lapply(
          idx_list,
          function(idx) balanced$A[idx, , drop = FALSE]
        ),
        dvec = balanced$d,
        param = param
      ))
    }

    row_blocks <- vector("list", L)
    k3_by_block <- vector("list", L)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      row_blocks[[i]] <- A[idx, , drop = FALSE]
      k3_by_block[[i]] <- .extract_K3_slices(
        K3fun = K3op_list[[i]],
        tvec = tvec[idx],
        param = param,
        block_dim = dims[i],
        basis = basis_list[[i]]
      )
    }

    .k3_slices_abcabc_from_factored_Q(
      k3_by_block = k3_by_block,
      row_blocks = row_blocks,
      dvec = dvec,
      param = param
    )
  }
  ineq_constraint <- function(tvec, param) {
    if (length(tvec) != total_dim) {
      stop(sprintf("`tvec` length mismatch in ineq_constraint: got %d, expected %d.",
                   length(tvec), total_dim))
    }

    pieces <- vector("list", L)
    for (i in seq_len(L)) {
      idx <- idx_list[[i]]
      pieces[[i]] <- ineq_list[[i]](tvec[idx], param)
    }
    total_size <- sum(lengths(pieces))

    out <- numeric(total_size)*param[1]
    if (total_size == 0) return(out)

    pos <- 1
    for (p in pieces) {
      lp <- length(p)
      if (lp > 0L) {
        out[pos:(pos + lp - 1)] <- p
        pos <- pos + lp
      }
    }
    out
  }


  # ---------------------------------------------------------------------------
  # Analytic t-hat (if all children have one)
  # ---------------------------------------------------------------------------
  analytic_tvec_hat <- NULL
  if (!is.null(analytic_hat_list)) {
    analytic_tvec_hat <- function(x, param) {
      if (length(x) != total_dim) {
        stop(sprintf("analytic_tvec_hat: `x` length %d, expected %d.", length(x), total_dim))
      }
      out <- numeric(total_dim) * param[1]
      for (i in seq_len(L)) {
        idx <- idx_list[[i]]
        out[idx] <- analytic_hat_list[[i]](x[idx], param)
      }
      out
    }
  }

  # call_history => for debugging
  hist_pieces <- vapply(
    cgf_list,
    function(cg) paste(cg$call_history, collapse = " -> "),
    character(1)
  )
  combined_history <- paste0("[", paste(hist_pieces, collapse = ", "), "]")
  op_name <- c(combined_history, "concatenationCGF")


  # simulation (only if all components can simulate)
  rsim <- NULL
  if (all(vapply(cgf_list, function(cg) isTRUE(cg$has_rsim), logical(1)))) {
    rsim_list <- lapply(cgf_list, function(cg) cg$rsim)
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {

      if (vector_length != total_dim) {
        stop(
          "concatenationCGF$rsim: 'vector_length' must equal sum(component_dims) = ",
          total_dim, ". Got ", vector_length, ".",
          call. = FALSE
        )
      }

      out <- matrix(0, nrow = total_dim, ncol = n)

      for (i in seq_len(L)) {
        idx <- idx_list[[i]]
        t_sub <- if (is.null(tvec)) NULL else tvec[idx]

        out[idx, ] <- rsim_list[[i]](
          n = n,
          vector_length = dims[i],
          parameter_vector = parameter_vector,
          tvec = t_sub,
          flatten = FALSE,
          ...
        )
      }

      out
    }
  }

  # Build args list (names match createCGF parameters exactly)
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
    op_name = op_name
  )

  extra_args <- list(...)
  extra_args <- extra_args[!vapply(extra_args, is.null, logical(1))]
  extra_names <- names(extra_args)

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
      "concatenationCGF cannot override generated method(s) through ",
      "'...': ", paste(conflicting_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (any(contraction_names %in% extra_names)) cgf_args$func_T <- NULL

  do.call(createCGF, modifyList(cgf_args, extra_args))
}






#' @title Concatenation of CGF objects
#'
#' @description
#' Constructs a `CGF` object for a concatenated random vector whose components are
#' independent and each specified by its own `CGF` object.
#'
#' If \code{cgf_list = list(C1, ..., CL)} and \code{component_dims = (d1,...,dL)},
#' then \code{tvec} is partitioned as
#' \deqn{tvec = (t^{(1)}, ..., t^{(L)}) \quad\text{with}\quad length(t^{(i)}) = d_i,}
#' and the resulting CGF satisfies
#' \deqn{K(tvec;\theta) = \sum_{i=1}^L K_i(t^{(i)};\theta).}
#'
#'
#' @param cgf_list A non-empty list of `CGF` objects.
#' @param component_dims An integer vector of length \code{length(cgf_list)} giving the
#'   dimension of each component. If a single integer is supplied, it is recycled.
#' @param iidReps Either \code{"any"} (default) or a positive integer. See Details.
#' @param ... Additional arguments passed to \code{createCGF()} (advanced use).
#'
#' @seealso \code{\link{iidReplicatesCGF}}, \code{\link{sumOfIndependentCGF}}
#'
#' @return A `CGF` object.
#'
#' @examples
#' ## univariate + univariate (Poisson + Gamma)
#' # Use adaptors so both components share a single global parameter vector:
#' #   theta = c(lambda, shape, rate)
#' cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
#' cg_gam  <- GammaModelCGF(shape = adaptor(indices = 2),
#'                          rate  = adaptor(indices = 3),
#'                          iidReps = 1)
#'
#' # component_dims can be a scalar (recycled) or a vector:
#' cg <- concatenationCGF(list(cg_pois, cg_gam),
#'                        component_dims = 1L,
#'                        iidReps = "any")
#'
#' theta <- c(lambda = 2, shape = 5, rate = 3)
#'
#' # Two iid blocks; each block is (Poisson, Gamma) so block_size = 2.
#' tvec <- c(0.05, 0.10,
#'          -0.02, 0.20)
#'
#' # K is additive across blocks and across components:
#' K_manual <- cg_pois$K(tvec[1], theta) + cg_gam$K(tvec[2], theta) +
#'             cg_pois$K(tvec[3], theta) + cg_gam$K(tvec[4], theta)
#' stopifnot(all.equal(cg$K(tvec, theta), K_manual))
#'
#' # K1 is blockwise concatenation:
#' K1_manual <- c(cg_pois$K1(tvec[1], theta), cg_gam$K1(tvec[2], theta),
#'                cg_pois$K1(tvec[3], theta), cg_gam$K1(tvec[4], theta))
#' stopifnot(all.equal(as.numeric(cg$K1(tvec, theta)),
#'                     as.numeric(K1_manual)))
#'
#' # Analytic t-hat is available if all children have it:
#' if (isTRUE(cg$has_analytic_tvec_hat)) {
#'   x <- c(12, 5, 8, 7)  # observed data (same layout as tvec)
#'   t_hat <- cg$analytic_tvec_hat(x, theta)
#'   t_hat_manual <- c(cg_pois$analytic_tvec_hat(x[1], theta),
#'                     cg_gam$analytic_tvec_hat(x[2], theta),
#'                     cg_pois$analytic_tvec_hat(x[3], theta),
#'                     cg_gam$analytic_tvec_hat(x[4], theta))
#'   stopifnot(all.equal(as.numeric(t_hat), as.numeric(t_hat_manual)))
#' }
#'
#' ## Multivariate + univariate (Multinomial + Poisson)
#' # Multinomial dimension d = 3; concatenate with a Poisson scalar -> block_size = 4.
#' cg_mult <- MultinomialModelCGF(
#'   n        = adaptor(fixed_param = 10),
#'   prob_vec = adaptor(indices = 2:4),
#'   iidReps  = 1
#' )
#' cg_pois2 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
#'
#' cg2 <- concatenationCGF(list(cg_mult, cg_pois2),
#'                         component_dims = c(3L, 1L),
#'                         iidReps = "any")
#'
#' theta2 <- c(lambda = 4, p1 = 0.2, p2 = 0.3, p3 = 0.5)
#'
#' # Two iid blocks, each of length 4:
#' tvec2 <- c(0.02, -0.01, 0.03,  0.10,
#'           -0.02,  0.00, 0.01, -0.05)
#'
#' # Shape checks:
#' stopifnot(length(cg2$K1(tvec2, theta2)) == length(tvec2))
#' K2_full <- cg2$K2(tvec2, theta2)
#' stopifnot(all(dim(K2_full) == c(length(tvec2), length(tvec2))))
#'
#'
#' ## find.saddlepoint.MLE on a concatenation model
#' \donttest{
#' set.seed(123)
#' B <- 30
#' lambda_true <- c(5, 15)
#' Y <- rbind(rpois(B, lambda_true[1]),
#'            rpois(B, lambda_true[2]))
#'
#' cg_1 <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
#' cg_2 <- PoissonModelCGF(lambda = adaptor(indices = 2), iidReps = 1)
#' cg_pp <- concatenationCGF(list(cg_1, cg_2),
#'                           component_dims = 1L,
#'                           iidReps = "any")
#'
#' fit <- find.saddlepoint.MLE(
#'   observed.data  = Y, # 2 x B matrix -> columns are iid blocks
#'   cgf            = cg_pp,
#'   starting.theta = c(0.1, 0.9),
#'   lb.theta       = c(0.001, 0.1),
#'   method         = "two_step"
#' )
#'
#' stopifnot(max(abs(fit$MLEs.theta - rowMeans(Y))) < 1e-5)
#' }
#' @export
concatenationCGF <- function(cgf_list,
                             component_dims = 1L,
                             iidReps = "any",
                             ...) {
  .check_iidReps(iidReps)

  if (!all(component_dims == as.integer(component_dims)) || any(component_dims <= 0)) {
    stop("'component_dims' must be positive integers.")
  } else if (length(component_dims) == 1) {
    component_dims <- rep.int(as.integer(component_dims), length(cgf_list))
  } else if (length(component_dims) != length(cgf_list)) {
    stop("'component_dims' must have the same length as 'cgf_list'.")
  }

  base_cgf <- .concatenationCGF_internal(
    cgf_list = cgf_list,
    component_dims = component_dims,
    ...
  )

  # For concatenation, the natural block size is fixed and known:
  # one "observation block" is one full concatenated vector.
  block_size <- sum(component_dims)

  iidReplicatesCGF(cgf = base_cgf, iidReps = iidReps, block_size = block_size)
}
