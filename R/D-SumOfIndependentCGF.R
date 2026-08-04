# ------------------------------------------------------------
#  R/SumOfIndependentCGF.R
#
#   We build a CGF for the sum of *independent* random vectors:
#      Y = Y^(1) + ... + Y^(L)
#    where each summand has its own CGF and all summands share the same t-vector.
#
#    - Cache method closures from each CGF in cgf_list to reduce R6 "$" dispatch.
#    - Override the *private* factored operator methods
#        K4operatorAABB_factored
#        K3K3operatorAABBCC_factored
#        K3K3operatorABCABC_factored
#      with symmetry-exploiting implementations for the common repeated-Q case.
#      These are used by the base CGF default correction term func_T().
# ------------------------------------------------------------

.sumOfIndependentCGF_internal <- function(cgf_list, ...) {

  # Cache methods
  K_list       <- lapply(cgf_list, function(cg) cg$K)
  K1_list      <- lapply(cgf_list, function(cg) cg$K1)
  K2_list      <- lapply(cgf_list, function(cg) cg$K2)
  K3op_list    <- lapply(cgf_list, function(cg) cg$K3operator)
  K4op_list    <- lapply(cgf_list, function(cg) cg$K4operator)

  K2op_list       <- lapply(cgf_list, function(cg) cg$K2operator)
  K2opAK2AT_list  <- lapply(cgf_list, function(cg) cg$K2operatorAK2AT)
  K2_factor_list  <- lapply(cgf_list, .K2_factor_method)
  K4AABB_list     <- lapply(cgf_list, function(cg) cg$K4operatorAABB)
  K3K3AABBCC_list <- lapply(cgf_list, function(cg) cg$K3K3operatorAABBCC)
  K3K3ABCABC_list <- lapply(cgf_list, function(cg) cg$K3K3operatorABCABC)
  K4AABB_factored_list <- lapply(
    cgf_list, function(cg) cg$.private_api$K4operatorAABB_factored
  )
  K3K3AABBCC_factored_list <- lapply(
    cgf_list, function(cg) cg$.private_api$K3K3operatorAABBCC_factored
  )
  K3K3ABCABC_factored_list <- lapply(
    cgf_list, function(cg) cg$.private_api$K3K3operatorABCABC_factored
  )
  K4AABB_delegate_safe <- vapply(
    K4AABB_factored_list, .factored_delegate_is_safe, logical(1)
  )
  K3K3AABBCC_delegate_safe <- vapply(
    K3K3AABBCC_factored_list, .factored_delegate_is_safe, logical(1)
  )
  K3K3ABCABC_delegate_safe <- vapply(
    K3K3ABCABC_factored_list, .factored_delegate_is_safe, logical(1)
  )

  func_T <- NULL
  if (length(cgf_list) == 1L) {
    child_func_T <- cgf_list[[1L]]$.private_api$func_T
    # Do not pass the child's bound R6 method directly to createCGF(): doing so
    # rebinds its `self` to the sum.  This forwarding closure retains the
    # child's authoritative method and environment.
    func_T <- local({
      child_T <- child_func_T
      function(tvec, param) child_T(tvec, param)
    })
  }

  validate_factored <- function(tvec, A, d, where) {
    A_dim <- dim(A)
    if (length(A_dim) != 2L || A_dim[1L] != length(tvec)) {
      stop(where, ": A must have nrow(A) == length(tvec).", call. = FALSE)
    }
    if (A_dim[2L] != length(d)) {
      stop(where, ": ncol(A) must equal length(d).", call. = FALSE)
    }
    length(d)
  }

  tilting_list <- lapply(cgf_list, function(cg) cg$.private_api$tilting_exponent)

  ineq_list <- lapply(cgf_list, function(cg) cg$ineq_constraint)




  K <- function(tvec, param) {
    total <- .ad_zero_scalar(param)
    for (f in K_list) total <- total + f(tvec, param)
    total
  }

  K1 <- function(tvec, param) {
    out <- numeric(length(tvec)) * .ad_type_scale(param)
    for (f in K1_list) out <- out + f(tvec, param)
    out
  }

  K2 <- function(tvec, param) {
    d <- length(tvec)
    accum <- matrix(0, nrow = d, ncol = d) * .ad_type_scale(param)
    for (f in K2_list) accum <- accum + f(tvec, param)
    accum
  }

  # Factored contractions reuse the state and certify all child sums with one
  # terminal atomic.  Public primitive calls still certify immediately; this
  # avoids an R atomic crossing inside every rank/coordinate loop.
  K3operator_sum_state <- function(tvec, param, v1, v2, v3) {
    total <- .ad_compensated_state(param)
    for (f in K3op_list) {
      total <- .ad_compensated_add(
        total,
        f(tvec, param, v1, v2, v3)
      )
    }
    total
  }

  K3operator <- function(tvec, param, v1, v2, v3) {
    if (length(K3op_list) == 1L) {
      return(K3op_list[[1L]](tvec, param, v1, v2, v3))
    }

    total <- K3operator_sum_state(tvec, param, v1, v2, v3)
    .ad_guard_compensated_sums(.ad_compensated_value(total), list(total))
  }

  K4operator_sum_state <- function(tvec, param, v1, v2, v3, v4) {
    total <- .ad_compensated_state(param)
    for (f in K4op_list) {
      total <- .ad_compensated_add(
        total,
        f(tvec, param, v1, v2, v3, v4)
      )
    }
    total
  }

  K4operator <- function(tvec, param, v1, v2, v3, v4) {
    if (length(K4op_list) == 1L) {
      return(K4op_list[[1L]](tvec, param, v1, v2, v3, v4))
    }

    total <- K4operator_sum_state(tvec, param, v1, v2, v3, v4)
    .ad_guard_compensated_sums(.ad_compensated_value(total), list(total))
  }

  tilting_exponent <- function(tvec, param) {
    total <- .ad_zero_scalar(param)
    for (f in tilting_list) total <- total + f(tvec, param)
    total
  }

  K2operator <- function(tvec, param, x, y) {
    total <- .ad_zero_scalar(param)
    for (f in K2op_list) total <- total + f(tvec, param, x, y)
    total
  }

  K2operatorAK2AT <- function(tvec, param, B) {
    r <- nrow(B)
    accum <- matrix(0, nrow = r, ncol = r) * .ad_type_scale(param)
    for (f in K2opAK2AT_list) accum <- accum + f(tvec, param, B)
    accum
  }

  K2_factor <- NULL
  if (any(vapply(K2_factor_list, is.function, logical(1)))) {
    K2_factor <- function(tvec, param, B) {
      terms <- list()
      for (i in seq_along(K2_factor_list)) {
        child_terms <- if (is.function(K2_factor_list[[i]])) {
          K2_factor_list[[i]](tvec, param, B)
        } else {
          .K2_dense_term(K2opAK2AT_list[[i]](tvec, param, B))
        }
        terms[[length(terms) + 1L]] <- child_terms
      }
      unlist(terms, recursive = FALSE)
    }
  }

  K4operatorAABB <- function(tvec, param, Q) {
    if (length(K4AABB_list) == 1L) {
      return(K4AABB_list[[1L]](tvec, param, Q))
    }

    total <- .ad_compensated_state(param)
    for (f in K4AABB_list) {
      total <- .ad_compensated_add(total, f(tvec, param, Q))
    }
    .ad_guard_compensated_sums(.ad_compensated_value(total), list(total))
  }

  K3K3operatorAABBCC <- if (length(cgf_list) == 1L) {
    local({
      child_method <- K3K3AABBCC_list[[1L]]
      function(tvec, param, Q) child_method(tvec, param, Q)
    })
  } else {
    NULL
  }

  K3K3operatorABCABC <- if (length(cgf_list) == 1L) {
    local({
      child_method <- K3K3ABCABC_list[[1L]]
      function(tvec, param, Q) child_method(tvec, param, Q)
    })
  } else {
    NULL
  }

  ineq_constraint <- function(tvec, param) {
    pieces <- lapply(ineq_list, function(f) f(tvec, param))
    total_size <- sum(lengths(pieces))

    out <- numeric(total_size) * .ad_type_scale(param)
    if (total_size == 0L) return(out)

    idx <- 1L
    for (p in pieces) {
      lp <- length(p)
      if (lp > 0) {
        out[idx:(idx + lp - 1L)] <- p
        idx <- idx + lp
      }
    }
    out
  }




  K4operatorAABB_factored <- function(tvec, param, A, d) {

    r <- validate_factored(
      tvec, A, d, "sumOfIndependent K4operatorAABB_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    if (length(cgf_list) == 1L) {
      return(K4AABB_factored_list[[1L]](tvec, param, A, d))
    }

    balanced <- .balance_factored_Q(
      A, d, "sumOfIndependent K4operatorAABB_factored"
    )
    A <- balanced$A
    d <- balanced$d

    Acols <- lapply(seq_len(r), function(i) A[, i])

    res <- .ad_compensated_state(param)
    diagnostics <- list()
    for (i in seq_len(r)) {
      ai <- Acols[[i]]
      di <- d[i]
      for (j in i:r) {
        aj <- Acols[[j]]
        mult <- if (i == j) 1 else 2
        child_sum <- K4operator_sum_state(
          tvec, param, ai, ai, aj, aj
        )
        diagnostics[[length(diagnostics) + 1L]] <- child_sum
        res <- .ad_compensated_add(
          res,
          mult * (di * d[j]) * .ad_compensated_value(child_sum)
        )
      }
    }
    .ad_guard_compensated_sums(
      .ad_compensated_value(res), c(diagnostics, list(res))
    )
  }


  K3K3operatorAABBCC_factored <- function(tvec, param, A, d) {

    r <- validate_factored(
      tvec, A, d, "sumOfIndependent K3K3operatorAABBCC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    if (length(cgf_list) == 1L && (
      K3K3AABBCC_delegate_safe[[1L]] ||
        .use_direct_factored_rank(length(tvec), r, "AABBCC")
    )) {
      return(K3K3AABBCC_factored_list[[1L]](
        tvec, param, A, d
      ))
    }

    if (length(cgf_list) == 1L) {
      balanced <- .balance_factored_Q(
        A, d, "sumOfIndependent K3K3operatorAABBCC_factored fallback"
      )
      A <- balanced$A
      d <- balanced$d
      slices <- .extract_K3_slices(
        K3operator, tvec, param, length(tvec), diag(1, length(tvec))
      )
      Q <- .factor_block_matrix(A, d, A)
      u <- .k3_slices_to_aabbcc_vector(slices, Q, param)
      z <- as.vector(crossprod(A, u))
      return(sum(d * z * z))
    }

    balanced <- .balance_factored_Q(
      A, d, "sumOfIndependent K3K3operatorAABBCC_factored"
    )
    A <- balanced$A
    d <- balanced$d

    Acols <- lapply(seq_len(r), function(i) A[, i])

    res <- .ad_compensated_state(param)
    diagnostics <- list()
    for (j in seq_len(r)) {
      aj <- Acols[[j]]
      g <- .ad_compensated_state(param)
      for (i in seq_len(r)) {
        ai <- Acols[[i]]
        child_sum <- K3operator_sum_state(tvec, param, ai, ai, aj)
        diagnostics[[length(diagnostics) + 1L]] <- child_sum
        g <- .ad_compensated_add(
          g,
          d[i] * .ad_compensated_value(child_sum)
        )
      }
      diagnostics[[length(diagnostics) + 1L]] <- g
      g_value <- .ad_compensated_value(g)
      res <- .ad_compensated_add(res, (d[j] * g_value) * g_value)
    }
    .ad_guard_compensated_sums(
      .ad_compensated_value(res), c(diagnostics, list(res))
    )
  }

  K3K3operatorABCABC_factored <- function(tvec, param, A, d) {

    # Symmetry exact path (i<=j<=k with multiplicities 1/3/6)
    r <- validate_factored(
      tvec, A, d, "sumOfIndependent K3K3operatorABCABC_factored"
    )
    if (r == 0L) return(.ad_zero_scalar(param))

    if (length(cgf_list) == 1L && (
      K3K3ABCABC_delegate_safe[[1L]] ||
        .use_direct_factored_rank(length(tvec), r, "ABCABC")
    )) {
      return(K3K3ABCABC_factored_list[[1L]](
        tvec, param, A, d
      ))
    }

    if (length(cgf_list) == 1L) {
      slices <- .extract_K3_slices(
        K3operator, tvec, param, length(tvec), diag(1, length(tvec))
      )
      return(.k3_slices_abcabc_from_factored_Q(
        list(slices), list(A), d, param
      ))
    }

    balanced <- .balance_factored_Q(
      A, d, "sumOfIndependent K3K3operatorABCABC_factored"
    )
    A <- balanced$A
    d <- balanced$d

    Acols <- lapply(seq_len(r), function(i) A[, i])

    res <- .ad_compensated_state(param)
    diagnostics <- list()
    for (i in seq_len(r)) {
      ai <- Acols[[i]]
      di <- d[i]
      for (j in i:r) {
        aj <- Acols[[j]]
        dij <- di * d[j]
        for (k in j:r) {
          ak <- Acols[[k]]
          child_sum <- K3operator_sum_state(tvec, param, ai, aj, ak)
          diagnostics[[length(diagnostics) + 1L]] <- child_sum
          val <- .ad_compensated_value(child_sum)

          mult <- if (i == j && j == k) {
            1
          } else if (i == j || j == k) {
            3
          } else {
            6
          }

          res <- .ad_compensated_add(
            res,
            mult * (dij * d[k]) * (val * val)
          )
        }
      }
    }
    .ad_guard_compensated_sums(
      .ad_compensated_value(res), c(diagnostics, list(res))
    )
  }
  # A singleton K4 delegate is bounded exactly when its child is.  Multi-child
  # K4 rank loops are deliberately marked unsafe so an outer consumer can use
  # the public coordinate contraction instead.
  K4operatorAABB_factored <- .factored_delegate_mark(
    K4operatorAABB_factored,
    length(cgf_list) == 1L && K4AABB_delegate_safe[[1L]]
  )

  # A singleton either delegates a safe child or uses the bounded coordinate
  # fallback above. Multi-child rank formulas do not select against coordinate
  # work and therefore remain unsafe to delegate through another singleton.
  K3K3operatorAABBCC_factored <- .factored_delegate_mark(
    K3K3operatorAABBCC_factored, length(cgf_list) == 1L
  )
  K3K3operatorABCABC_factored <- .factored_delegate_mark(
    K3K3operatorABCABC_factored, length(cgf_list) == 1L
  )

  # -------------------------------------------------------------------------
  # Call history label (must collapse each child's call_history first)
  # -------------------------------------------------------------------------
  hist_pieces <- vapply(
    cgf_list,
    function(cg) paste(cg$call_history, collapse = " -> "),
    character(1)
  )
  combined_history <- paste0("[", paste(hist_pieces, collapse = ", "), "]")
  op_name_vec <- c(combined_history, "sumOfIndependentCGF")

  # simulation (only if all summands can simulate)
  rsim <- NULL
  if (all(vapply(cgf_list, function(cg) isTRUE(cg$has_rsim), logical(1)))) {
    rsim_list <- lapply(cgf_list, function(cg) cg$rsim)
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      out <- rsim_list[[1]](
        n = n,
        vector_length = vector_length,
        parameter_vector = parameter_vector,
        tvec = tvec,
        flatten = FALSE,
        ...
      )

      if (length(rsim_list) > 1L) {
        for (j in 2:length(rsim_list)) {
          out <- out + rsim_list[[j]](
            n = n,
            vector_length = vector_length,
            parameter_vector = parameter_vector,
            tvec = tvec,
            flatten = FALSE,
            ...
          )
        }
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
    ineq_constraint = ineq_constraint,
    K2operator = K2operator,
    K2operatorAK2AT = K2operatorAK2AT,
    K2_factor = K2_factor,
    K2_factor_terminal = if (!is.null(K2_factor)) function() TRUE else NULL,
    func_T = func_T,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    rsim = rsim,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    op_name = op_name_vec
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
      "sumOfIndependentCGF cannot override generated method(s) through ",
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



#' @title CGF Object for the sum of independent random variables
#'
#' @description
#' Constructs a new CGF object representing the sum of independent random vectors:
#' \eqn{Y = Y^{(1)} + \cdots + Y^{(L)}} where each summand has its own CGF and the summands are independent.
#'
#' @details
#' **Replication (iidReps / block_size):**
#' \itemize{
#'   \item If both `iidReps` and `block_size` are `NULL`, no replication is applied.
#'   \item If `block_size` is provided but `iidReps` is `NULL`, we set `iidReps = "any"` and
#'         infer the number of blocks from `length(tvec) / block_size` at evaluation time.
#'   \item If `iidReps = "any"`, then `block_size` must be provided.
#'   \item If `iidReps` is a positive integer, `block_size` may be `NULL`, though providing it is encouraged.
#' }
#'
#' Note: \code{iidReps}/\code{block_size} describe i.i.d. replication of the sum \eqn{Y}, not the length of \code{cgf_list}.
#'
#' @param cgf_list A non-empty list of CGF objects (each inherits from class \code{"CGF"}).
#' @param iidReps Optional. Either `NULL`, \code{"any"}, or a positive integer.
#' @param block_size Either \code{NULL} or a positive integer describing the block size for iid replication.
#' @param ... Additional named arguments passed to \code{\link{createCGF}} (rare).
#'
#' @examples
#' ## -----------------------------
#' ## Sum of independent Poisson variables
#' ## -----------------------------
#' ## Here we build Y = Y1 + Y2 with Y1 ~ Pois(lambda1), Y2 ~ Pois(lambda2),
#' ## independent. The mean is lambda1 + lambda2, i.e. K1(0) should equal that.
#' lambda1 <- 2
#' lambda2 <- 5
#'
#' cg1 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda1), iidReps = 1)
#' cg2 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda2), iidReps = 1)
#'
#' ## one observation (block_size = 1)
#' cg_sum <- sumOfIndependentCGF(list(cg1, cg2), iidReps = 1, block_size = 1)
#'
#' theta_dummy <- 0
#' cg_sum$K1(tvec = 0, parameter_vector = theta_dummy)   # should be 7
#'
#'
#' ## ------------------------------------------------------------
#' ## Sum of independent Gammas with different rates
#' ## (demonstrates concatenated inequality constraints)
#' ## ------------------------------------------------------------
#' ## If X ~ Gamma(shape=a, rate=b), its CGF is only valid for tvec < b.
#' ## For Y = X1 + X2 with rates b1 and b2, the domain is t < min(b1, b2),
#' ## and the implementation is such that the constraint vector is concatenated:
#' ##   g(t,theta) = c(t - b1, t - b2)  <= 0
#' ##
#' \donttest{
#' set.seed(123)
#' B <- 50
#'
#' ## True shapes (fixed) and rates (unknown; theta = (b1, b2))
#' a1_true <- 10
#' a2_true <- 3
#' b1_true <- 7
#' b2_true <- 1
#' theta_true <- c(b1_true, b2_true)
#'
#' ## Two GammaModelCGFs, each uses one component of theta as its rate
#' cg_g1 <- GammaModelCGF(
#'   shape  = adaptor(fixed_param = a1_true),
#'   rate   = function(th) th[1],
#'   iidReps = 1
#' )
#' cg_g2 <- GammaModelCGF(
#'   shape  = adaptor(fixed_param = a2_true),
#'   rate   = function(th) th[2],
#'   iidReps = 1
#' )
#'
#' ## Sum for ONE observation
#' cg_sum_one <- sumOfIndependentCGF(list(cg_g1, cg_g2), iidReps = 1, block_size = 1)
#'
#' ## Inequality constraints are concatenated across summands:
#' ## For scalar t: g(t,theta) = c(t-b1, t-b2) must be <= 0.
#' cg_sum_one$ineq_constraint(tvec = 0.9, param = theta_true)  # ~ c(-6.1, -0.1) (feasible)
#' cg_sum_one$ineq_constraint(tvec = 1.2, param = theta_true)  # ~ c(-5.8, +0.2) (violates 2nd)
#'
#' ## B i.i.d. replicates of the sum
#' cg_sum_B <- sumOfIndependentCGF(list(cg_g1, cg_g2), iidReps = B, block_size = 1)
#'
#' ## Simulate data y_i = x1_i + x2_i
#' y <- rgamma(B, shape = a1_true, rate = b1_true) +
#'      rgamma(B, shape = a2_true, rate = b2_true)
#'
#' ## NOTE / recommendation:
#' ## When cgf$ineq_constraint is non-empty (domain-constrained CGFs),
#' ## method="constrained" is currently the most robust choice because it enforces
#' ## the CGF domain constraints directly during optimisation. The "two_step" method
#' ## WILL be slower for constrained CGFs in the current implementation.
#'
#' fit_const <- find.saddlepoint.MLE(
#'   observed.data  = y,
#'   cgf            = cg_sum_B,
#'   starting.theta = c(1.2, 0.5),
#'   lb.theta       = c(1e-4, 2e-5),
#'   method         = "constrained"
#' )
#'
#' fit_two <- find.saddlepoint.MLE(
#'   observed.data  = y,
#'   cgf            = cg_sum_B,
#'   starting.theta = c(1.2, 0.5),
#'   lb.theta       = c(1e-4, 2e-5),
#'   method         = "two_step"
#' )
#'
#' ## Optional quick checks: saddlepoint residual and feasibility
#' max(abs(cg_sum_B$K1(fit_const$MLEs.tvec, fit_const$MLEs.theta) - y))
#' max(cg_sum_B$ineq_constraint(fit_const$MLEs.tvec, fit_const$MLEs.theta))  # should be <= 0
#'
#' cat("true theta:", theta_true, "\n")
#' cat("constrained:", round(fit_const$MLEs.theta, 4), "\n")
#' cat("two_step   :", round(fit_two$MLEs.theta, 4), "\n")
#' }
#'
#'
#' @return A `CGF` object.
#' @export
sumOfIndependentCGF <- function(cgf_list,
                                iidReps = NULL,
                                block_size = NULL,
                                ...) {
  if (!is.list(cgf_list) || length(cgf_list) == 0) stop("'cgf_list' must be a non-empty list of CGF objects.")
  if (any(vapply(cgf_list, function(x) !inherits(x, "CGF"), FALSE)) ) stop("Every element of 'cgf_list' must be of class 'CGF'.")

  base_cgf <- .sumOfIndependentCGF_internal(cgf_list, ...)

  if (is.null(block_size) && is.null(iidReps)) return(base_cgf)
  if (is.null(iidReps)) iidReps <- "any"
  .check_iidReps(iidReps)
  if (is.numeric(iidReps) && iidReps == 1L && is.null(block_size)) return(base_cgf)
  if (identical(iidReps, "any") && is.null(block_size)) {
    stop("sumOfIndependentCGF(): iidReps='any' requires a non-NULL 'block_size'.")
  }

  iidReplicatesCGF(cgf = base_cgf, iidReps = iidReps, block_size = block_size)
}
