# ------------------------------------------------------------
# R/D-ShiftedCGF.R
#
# deterministic shift of a CGF
#   Y = X + b(theta)      (b may be constant or theta-dependent)
#
# Then:
#   K_Y(t;theta)  = K_X(t;theta) + <t, b(theta)>
#   K1_Y(t;theta) = K1_X(t;theta) + b(theta)
#   K2_Y, K3_Y, K4_Y are unchanged from X
#
# Domain/ineq constraints in t are unchanged.
# tilting_exponent, neg_ll, func_T are invariant under shifts (for fixed t).
# Analytic t-hat (if available) transforms as:
#   t_hat_Y(y;theta) = t_hat_X(y - b(theta); theta)
# ------------------------------------------------------------

#' @noRd
.expand_to_length_strict <- function(x, n, what = "vector") {
  x <- as.vector(x)
  # x <- x[]
  lx <- length(x)
  if (lx == 0L) stop(sprintf("%s must have positive length.", what))
  if (lx == n) return(x)
  if (lx == 1L) return(x[rep.int(1L, n)])
  if (n %% lx != 0L) {
    stop(sprintf(
      "%s length (%d) must be 1, length(tvec) (%d), or a divisor of length(tvec).",
      what, lx, n
    ))
  }
  x[rep.int(seq_len(lx), times = n %/% lx)]
}

#' @noRd
.shiftedCGF_internal <- function(base_cgf, shift_fn, ...) {
  stopifnot(inherits(base_cgf, "CGF"))
  stopifnot(is.function(shift_fn))

  # Cache base methods (avoid repeated R6 dispatch)
  K0   <- base_cgf$K
  K10  <- base_cgf$K1
  K20  <- base_cgf$K2
  K30  <- base_cgf$K3operator
  K40  <- base_cgf$K4operator

  K2op0      <- base_cgf$K2operator
  K2opAK2AT0 <- base_cgf$K2operatorAK2AT
  K2solve0   <- base_cgf$K2_solve
  logdet0    <- base_cgf$logdetK2

  K4AABB0    <- base_cgf$K4operatorAABB
  K3K3A0     <- base_cgf$K3K3operatorAABBCC
  K3K3B0     <- base_cgf$K3K3operatorABCABC

  ineq0 <- base_cgf$ineq_constraint

  # Private methods that are *invariant* under shifts (safe to reuse)
  tilt0  <- base_cgf$.get_private_method("tilting_exponent")
  negll0 <- base_cgf$.get_private_method("neg_ll")
  funcT0 <- base_cgf$.get_private_method("func_T")

  # Factored-operator private methods (also invariant under shifts)
  K4AABB_fact0 <- base_cgf$.get_private_method("K4operatorAABB_factored")
  K3K3A_fact0  <- base_cgf$.get_private_method("K3K3operatorAABBCC_factored")
  K3K3B_fact0  <- base_cgf$.get_private_method("K3K3operatorABCABC_factored")

  # Helper: b(theta) expanded to length(tvec)
  b_at <- function(theta, n) {
    b <- shift_fn(theta)
    .expand_to_length_strict(b, n, what = "shift")
  }

  Kfun <- function(tvec, param) {
    b <- b_at(param, length(tvec))
    K0(tvec, param) + sum(tvec * b)
  }

  K1fun <- function(tvec, param) {
    b <- b_at(param, length(tvec))
    K10(tvec, param) + b
  }

  # Unchanged derivatives/operators
  K2fun   <- function(tvec, param) K20(tvec, param)
  K3opfun <- function(tvec, param, v1, v2, v3) K30(tvec, param, v1, v2, v3)
  K4opfun <- function(tvec, param, v1, v2, v3, v4) K40(tvec, param, v1, v2, v3, v4)

  # Unchanged constraints in t
  ineqfun <- function(tvec, param) ineq0(tvec, param)

  # Unchanged “fast” methods
  K2opfun <- function(tvec, param, x, y) K2op0(tvec, param, x, y)
  K2opAK2ATfun <- function(tvec, param, B) K2opAK2AT0(tvec, param, B)

  K2solve_fun <- function(tvec, param, rhs) K2solve0(tvec, param, rhs)
  logdet_fun  <- function(tvec, param)      logdet0(tvec, param)

  K4AABB_fun   <- function(tvec, param, Q1, Q2)      K4AABB0(tvec, param, Q1, Q2)
  K3K3A_fun    <- function(tvec, param, Q1, Q2, Q3)  K3K3A0(tvec, param, Q1, Q2, Q3)
  K3K3B_fun    <- function(tvec, param, Q1, Q2, Q3)  K3K3B0(tvec, param, Q1, Q2, Q3)


  tilting_fun <- function(tvec, param) tilt0(tvec, param)
  negll_fun   <- function(tvec, param) negll0(tvec, param)
  funcT_fun   <- function(tvec, param) funcT0(tvec, param)

  #
  analytic_tvec_hat_func <- NULL
  if (isTRUE(base_cgf$has_analytic_tvec_hat)) {
    hat0 <- base_cgf$analytic_tvec_hat
    analytic_tvec_hat_func <- function(x, param) {
      b <- b_at(param, length(x))
      hat0(x - b, param)
    }
  }

  # Simulation: if X can simulate, Y = X + b(theta) can simulate by shifting draws.
  simulate_fun <- NULL
  if (isTRUE(base_cgf$has_simulate)) {
    simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      X_sim <- base_cgf$rsim(
        n = n,
        vector_length = vector_length,
        parameter_vector = parameter_vector,
        tvec = tvec,
        flatten = FALSE,
        ...
      )
      b <- b_at(parameter_vector, vector_length)
      X_sim + b
    }
  }

  #
  base_hist <- paste(base_cgf$call_history, collapse = " -> ")
  op_name_vec <- c(base_hist, "shiftedCGF")

  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K3operator = K3opfun,
    K4operator = K4opfun,

    # invariant pieces
    tilting_exponent = tilting_fun,
    neg_ll = negll_fun,
    func_T = funcT_fun,
    rsim = simulate_fun,

    ineq_constraint = ineqfun,
    analytic_tvec_hat = analytic_tvec_hat_func,

    # pass-through “fast” linear algebra / operators
    K2operator      = K2opfun,
    K2operatorAK2AT = K2opAK2ATfun,
    K2_solve        = K2solve_fun,
    logdetK2        = logdet_fun,

    K4operatorAABB      = K4AABB_fun,
    K3K3operatorAABBCC  = K3K3A_fun,
    K3K3operatorABCABC  = K3K3B_fun,

    # pass-through factored operators (important for speed in func_T pipelines)
    K4operatorAABB_factored     = K4AABB_fact0,
    K3K3operatorAABBCC_factored = K3K3A_fact0,
    K3K3operatorABCABC_factored = K3K3B_fact0,

    op_name = op_name_vec,
    ...
  )
}

#' Shifted CGF (deterministic translation)
#'
#' @description
#' Returns the CGF of \eqn{Y = X + b(\theta)}, where \eqn{X} has CGF \code{cgf} and
#' \eqn{b(\theta)} is a deterministic shift (constant or theta-dependent).
#'
#' This operation does not change the dimension of the random vector and therefore
#' does not impose any new i.i.d. replication semantics. If you want a replicated
#' interpretation of \code{tvec}, wrap the result with \code{\link{iidReplicatesCGF}}.
#'
#' @param cgf A \code{CGF} object.
#' @param shift A numeric scalar/vector, an \code{adaptor}, or a function \code{shift(theta)}
#'   returning the shift vector. If its length is 1, it is recycled; if its length
#'   divides \code{length(tvec)}, it is recycled blockwise.
#' @param ... Passed to \code{createCGF()} (advanced).
#'
#' @return A \code{CGF} object.
#'
#' @examples
#' ## Example 1: Poisson shifted by a constant (exact relationship)
#' lam <- 2
#' b   <- 3
#' cgS <- shiftedCGF(PoissonCGF, shift = b)
#' tt   <- 0.2
#' stopifnot(all.equal(cgS$K(tt, lam), PoissonCGF$K(tt, lam) + b*tt, tol = 1e-12))
#' stopifnot(all.equal(cgS$K1(tt, lam), PoissonCGF$K1(tt, lam) + b, tol = 1e-12))
#'
#' ## Example 2: Multivariate Normal with theta-dependent shift + MLE
#' \dontrun{
#' set.seed(1)
#' d <- 2
#' B <- 200
#'
#' Sigma <- matrix(c(1.0, 0.3,
#'                   0.3, 2.0), nrow = d, byrow = TRUE)
#' mu0 <- c(-1, 0.2)
#'
#' # Base MVN with fixed mean/covariance (use fixed_param to keep AD context)
#' cg_base <- MultivariateNormalModelCGF(
#'   mu    = adaptor(fixed_param = mu0),
#'   sigma = adaptor(fixed_param = Sigma),
#'   iidReps = "any"
#' )
#'
#' # Shift b(theta) in R^d
#' b_fun <- function(theta) c(theta[1], exp(theta[2]))
#' cg <- shiftedCGF(cg_base, shift = b_fun)
#'
#' # Simulate Y = X + b(theta_true)
#' theta_true <- c(0.7, log(2))     # b(theta_true) = (0.7, 2)
#' mu_true <- mu0 + b_fun(theta_true)
#' Y <- cg$rsim(n = B, vector_length = d, parameter_vector = theta_true)
#'
#' fit <- find.saddlepoint.MLE(
#'   observed.data  = Y,            # columns treated as i.i.d. blocks
#'   cgf            = cg,
#'   starting.theta = c(0, 0),
#'   method         = "two_step"
#' )
#' fit$MLEs.theta
#'
#' # Closed-form MLE (known Sigma): match the sample mean
#' m <- rowMeans(Y)
#' c(m[1] - mu0[1], log(m[2] - mu0[2]))
#' }
#'
#'
#' @export
shiftedCGF <- function(cgf,
                       shift,
                       ...) {
  stopifnot(inherits(cgf, "CGF"))

  dots <- list(...)
  if ("iidReps" %in% names(dots) || "block_size" %in% names(dots)) {
    stop("shiftedCGF() no longer accepts 'iidReps' or 'block_size'. ",
         "Wrap the result with iidReplicatesCGF(), e.g. ",
         "iidReplicatesCGF(shiftedCGF(cgf, shift = ...), iidReps = B, block_size = d).")
  }

  # Allow numeric shift directly
  if (is.numeric(shift) && !is.function(shift)) {
    shift <- adaptor(fixed_param = as.numeric(shift))
  }
  shift_fn <- validate_function_or_adaptor(shift)

  do.call(.shiftedCGF_internal, c(list(base_cgf = cgf, shift_fn = shift_fn), dots))
}
