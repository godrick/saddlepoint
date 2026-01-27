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
#      with symmetry-exploiting implementations for the common Q1=Q2=Q3 case.
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
  K4AABB_list     <- lapply(cgf_list, function(cg) cg$K4operatorAABB)

  tilting_list <- lapply(cgf_list, function(cg) cg$.private_api$tilting_exponent)

  ineq_list <- lapply(cgf_list, function(cg) cg$ineq_constraint)




  Kfun <- function(tvec, param) {
    total <- 0*param[1]
    for (f in K_list) total <- total + f(tvec, param)
    total
  }

  K1fun <- function(tvec, param) {
    out <- numeric(length(tvec)) * param[1]
    for (f in K1_list) out <- out + f(tvec, param)
    out
  }

  K2fun <- function(tvec, param) {
    d <- length(tvec)
    accum <- matrix(0, nrow = d, ncol = d) * param[1]
    # print(class(accum))
    for (f in K2_list) accum <- accum + f(tvec, param)
    accum
  }

  K3opfun <- function(tvec, param, v1, v2, v3) {
    total <- 0*param[1]
    for (f in K3op_list) total <- total + f(tvec, param, v1, v2, v3)
    total
  }

  K4opfun <- function(tvec, param, v1, v2, v3, v4) {
    total <- 0*param[1]
    for (f in K4op_list) total <- total + f(tvec, param, v1, v2, v3, v4)
    total
  }

  tiltingfun <- function(tvec, param) {
    total <- 0*param[1]
    for (f in tilting_list) total <- total + f(tvec, param)
    total
  }

  K2opfun <- function(tvec, param, x, y) {
    total <- 0*param[1]
    for (f in K2op_list) total <- total + f(tvec, param, x, y)
    total
  }

  K2opAK2ATfun <- function(tvec, param, B) {
    r <- nrow(B)
    accum <- matrix(0, nrow = r, ncol = r) * param[1]
    for (f in K2opAK2AT_list) accum <- accum + f(tvec, param, B)
    accum
  }

  K4AABBfun <- function(tvec, param, Q1, Q2) {
    total <- 0*param[1]
    for (f in K4AABB_list) total <- total + f(tvec, param, Q1, Q2)
    total
  }

  #

  ineqfun <- function(tvec, param) {
    pieces <- lapply(ineq_list, function(f) f(tvec, param))
    total_size <- sum(lengths(pieces))

    out <- numeric(total_size) * param[1]
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




  # -------------------------------------------------------------------------
  # Symmetric factored operators used by base func_T():
  #     K4operatorAABB_factored(t,theta,A,d,A,d)
  #     K3K3operatorAABBCC_factored(t,theta,A,d,A,d,A,d)
  #     K3K3operatorABCABC_factored(t,theta,A,d,A,d,A,d)
  # when Q1=Q2=Q3 (same A and d)
  # -------------------------------------------------------------------------

  K4operatorAABB_factored_sym <- function(tvec, param, A1, d1, A2, d2) {

    same_case <- identical(A1, A2) && identical(d1, d2)

    if (!same_case) {
      # generic fallback (matches base default)
      r1 <- length(d1); r2 <- length(d2)
      res <- 0 * param[1]
      for (m1 in seq_len(r1)) {
        for (m2 in seq_len(r2)) {
          res <- res + d1[m1] * d2[m2] * K4opfun(
            tvec, param, A1[, m1], A1[, m1], A2[, m2], A2[, m2]
          )
        }
      }
      return(res)
    }

    r <- length(d1)
    if (r == 0L) return(0*param[1])

    Acols <- lapply(seq_len(r), function(i) A1[, i])

    res <- 0*param[1]
    for (i in seq_len(r)) {
      ai <- Acols[[i]]
      di <- d1[i]
      for (j in i:r) {
        aj <- Acols[[j]]
        mult <- if (i == j) 1 else 2
        res <- res + mult * (di * d1[j]) * K4opfun(tvec, param, ai, ai, aj, aj)
      }
    }
    res
  }


  K3K3operatorAABBCC_factored_sym <- function(tvec, param, A1, d1, A2, d2, A3, d3) {

    same_case <- identical(A1, A2) && identical(A1, A3) &&
      identical(d1, d2) && identical(d1, d3)

    if (!same_case) {
      # generic fallback (base default)
      r1 <- length(d1); r2 <- length(d2); r3 <- length(d3)
      res <- 0*param[1]
      for (m2 in seq_len(r2)) {
        factor1 <- 0 * param[1]
        for (m1 in seq_len(r1)) {
          factor1 <- factor1 + d1[m1] * K3opfun(tvec, param, A1[, m1], A1[, m1], A2[, m2])
        }
        factor2 <- 0*param[1]
        for (m3 in seq_len(r3)) {
          factor2 <- factor2 + d3[m3] * K3opfun(tvec, param, A2[, m2], A3[, m3], A3[, m3])
        }
        res <- res + d2[m2] * factor1 * factor2
      }
      return(res)
    }

    # same-case exact simplification:
    #   res = sum_j d[j] * ( sum_i d[i] K3(a_i,a_i,a_j) )^2
    r <- length(d1)
    if (r == 0L) return(0*param[1])

    Acols <- lapply(seq_len(r), function(i) A1[, i])

    res <- 0*param[1]
    for (j in seq_len(r)) {
      aj <- Acols[[j]]
      g  <- 0 * param[1]
      for (i in seq_len(r)) {
        ai <- Acols[[i]]
        g <- g + d1[i] * K3opfun(tvec, param, ai, ai, aj)
      }
      res <- res + d1[j] * (g*g)
    }
    res
  }

  # ---- K3K3operatorABCABC_factored ----
  K3K3operatorABCABC_factored_sym <- function(tvec, param, A1, d1, A2, d2, A3, d3) {

    same_case <- identical(A1, A2) && identical(A1, A3) &&
      identical(d1, d2) && identical(d1, d3)

    if (!same_case) {
      # generic fallback (base default)
      r1 <- length(d1); r2 <- length(d2); r3 <- length(d3)
      res <- 0 * param[1]
      for (m1 in seq_len(r1)) {
        for (m2 in seq_len(r2)) {
          for (m3 in seq_len(r3)) {
            val <- K3opfun(tvec, param, A1[, m1], A2[, m2], A3[, m3])
            res <- res + d1[m1] * d2[m2] * d3[m3] * (val * val)
          }
        }
      }
      return(res)
    }

    # Symmetry exact path (i<=j<=k with multiplicities 1/3/6)
    r <- length(d1)
    if (r == 0L) return(0 * param[1])

    Acols <- lapply(seq_len(r), function(i) A1[, i])

    res <- 0 * param[1]
    for (i in seq_len(r)) {
      ai <- Acols[[i]]
      di <- d1[i]
      for (j in i:r) {
        aj <- Acols[[j]]
        dij <- di * d1[j]
        for (k in j:r) {
          ak <- Acols[[k]]
          val <- K3opfun(tvec, param, ai, aj, ak)

          mult <- if (i == j && j == k) {
            1
          } else if (i == j || j == k) {
            3
          } else {
            6
          }

          res <- res + mult * (dij * d1[k]) * (val * val)
        }
      }
    }
    res
  }

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
  simulate_fun <- NULL
  if (all(vapply(cgf_list, function(cg) isTRUE(cg$has_rsim), logical(1)))) {
    rsim_list <- lapply(cgf_list, function(cg) cg$rsim)
    simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
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


  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,

    tilting_exponent = tiltingfun,
    ineq_constraint  = ineqfun,

    K2operator       = K2opfun,
    K2operatorAK2AT  = K2opAK2ATfun,
    K4operatorAABB   = K4AABBfun,

    rsim = simulate_fun,

    # overrides for func_T path
    K4operatorAABB_factored     = K4operatorAABB_factored_sym,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored_sym,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored_sym,

    op_name = op_name_vec,
    ...
  )
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
