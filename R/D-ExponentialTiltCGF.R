# ------------------------------------------------------------------
#  R/ExponentialTiltCGF.R
#  Main function: ExponentialTiltCGF (exported)
#
#  -------
#  Build a CGF object for the exponentially tilted version of a base CGF:
#     dP_h(x) = exp(h^T x) dP(x) / E[exp(h^T X)]
#  If X has CGF K_X(t;theta), then under the exponential tilt with vector h:
#     K_{X,h}(t;theta) = K_X(t + h;theta) - K_X(h;theta).
#  --------------
#  - Same dimension as the base CGF (tvec dimension unchanged).
#  - Works for *any* CGF object (including ones produced by other operations).
#  - Domain / inequality constraints:
#       the tilted CGF requires BOTH:
#         base constraints at (t + h)  AND  base constraints at h.
#    So we concatenate: g_tilt(t) = c( g_base(t+h), g_base(h) ).
#
#  Replication / blocking
#  ----------------------
#  We follow the same conventions as sumOfIndependentCGF():
#    - The wrapper passes (iidReps, block_size) to iidReplicatesCGF().
#    - The tilt vector h(theta) is "recycled by blocks" automatically:
#        if iidReplicatesCGF is used, h(theta) is interpreted as a single-block tilt and repeated across blocks.
#        if length(h) == length(tvec), it is used as-is.
# ------------------------------------------------------------------

.exponentialTiltCGF_internal <- function(base_cgf, tilt_fn, ...) {
  stopifnot(inherits(base_cgf, "CGF"))
  stopifnot(is.function(tilt_fn))

  # Cache base methods to avoid repeated R6 dispatch
  K0   <- base_cgf$K
  K10  <- base_cgf$K1
  K20  <- base_cgf$K2
  K3o0 <- base_cgf$K3operator
  K4o0 <- base_cgf$K4operator

  K2op0      <- base_cgf$K2operator
  K2opAK2AT0 <- base_cgf$K2operatorAK2AT

  K2solve0 <- base_cgf$K2_solve
  logdet0  <- base_cgf$logdetK2

  K4AABB0     <- base_cgf$K4operatorAABB
  K3K3AABBCC0 <- base_cgf$K3K3operatorAABBCC
  K3K3ABCABC0 <- base_cgf$K3K3operatorABCABC

  ineq0 <- base_cgf$ineq_constraint

  # Factored private operators used by default func_T() implementations
  K4AABB_fact0     <- base_cgf$.private_api$K4operatorAABB_factored
  K3K3AABBCC_fact0 <- base_cgf$.private_api$K3K3operatorAABBCC_factored
  K3K3ABCABC_fact0 <- base_cgf$.private_api$K3K3operatorABCABC_factored

  # Analytic t-hat: if the base has one, we can shift it back by h
  has_analytic <- isTRUE(base_cgf$has_analytic_tvec_hat)
  hat0 <- if (has_analytic) base_cgf$analytic_tvec_hat else NULL

  # Helper: expand h(theta) to length(tvec)
  expand_h <- function(h, m) {
    if (length(h) == 0L) stop("ExponentialTiltCGF: tilt_fn(theta) returned length 0.")
    if (length(h) == m) return(h)
    if (m %% length(h) != 0L) {
      stop(
        "ExponentialTiltCGF: length(tvec) = ", m,
        " is not a multiple of length(h) = ", length(h), ".\n",
        "Provide a tilt vector h of length 1, block_size, or length(tvec)."
      )
    }
    h[(seq_len(m) - 1L) %% length(h) + 1L]
  }

  # ---- Core tilted cumulants ----
  K <- function(tvec, param) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    # type anchor via 0*param[1] pattern (RTMB friendliness)
    (K0(tvec + h, param) - K0(h, param)) + 0 * param[1]
  }

  K1 <- function(tvec, param) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K10(tvec + h, param) + 0 * param[1]
  }

  K2 <- function(tvec, param) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K20(tvec + h, param) + 0 * param[1]
  }

  K3operator <- function(tvec, param, v1, v2, v3) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K3o0(tvec + h, param, v1, v2, v3) + 0 * param[1]
  }

  K4operator <- function(tvec, param, v1, v2, v3, v4) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K4o0(tvec + h, param, v1, v2, v3, v4) + 0 * param[1]
  }

  # ---- Derived operators (shifted) ----
  K2operator <- function(tvec, param, x, y) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K2op0(tvec + h, param, x, y) + 0 * param[1]
  }

  K2operatorAK2AT <- function(tvec, param, Bmat) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K2opAK2AT0(tvec + h, param, Bmat) + 0 * param[1]
  }

  K2_solve <- function(tvec, param, rhs) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K2solve0(tvec + h, param, rhs) + 0 * param[1]
  }

  logdetK2 <- function(tvec, param) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    logdet0(tvec + h, param) + 0 * param[1]
  }

  K4operatorAABB <- function(tvec, param, Q) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K4AABB0(tvec + h, param, Q) + 0 * param[1]
  }

  K3K3operatorAABBCC <- function(tvec, param, Q) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K3K3AABBCC0(tvec + h, param, Q) + 0 * param[1]
  }

  K3K3operatorABCABC <- function(tvec, param, Q) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K3K3ABCABC0(tvec + h, param, Q) + 0 * param[1]
  }

  # Factored private operator wrappers (used by base func_T defaults)
  K4operatorAABB_factored <- function(tvec, param, A, d) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K4AABB_fact0(tvec + h, param, A, d) + 0 * param[1]
  }
  K3K3operatorAABBCC_factored <- function(tvec, param, A, d) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K3K3AABBCC_fact0(tvec + h, param, A, d) + 0 * param[1]
  }
  K3K3operatorABCABC_factored <- function(tvec, param, A, d) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)
    K3K3ABCABC_fact0(tvec + h, param, A, d) + 0 * param[1]
  }

  # Inequality constraints: g_base(t+h) and g_base(h)
  ineq_constraint <- function(tvec, param) {
    m <- length(tvec)
    h <- expand_h(tilt_fn(param), m)

    g1 <- ineq0(tvec + h, param)
    g0 <- ineq0(h, param)

    out <- numeric(length(g1) + length(g0)) * param[1]
    if (length(g1)) out[seq_along(g1)] <- g1
    if (length(g0)) out[length(g1) + seq_along(g0)] <- g0
    out
  }

  # Analytic t-hat shift: t_hat_tilt(x) = t_hat_base(x) - h
  analytic_tvec_hat <- NULL
  if (has_analytic) {
    analytic_tvec_hat <- function(x, param) {
      m <- length(x)
      h <- expand_h(tilt_fn(param), m)
      hat0(x, param) - h
    }
  }

  # This does NOT imply a general exact sampler from the tilted law.
  # rsim is only available if the base CGF can simulate under a tilt (tvec != NULL).
  rsim <- NULL
  if (isTRUE(base_cgf$has_rsim)) {
    rsim <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      h <- expand_h(tilt_fn(parameter_vector), vector_length)
      t_total <- if (is.null(tvec)) h else (h + tvec)

      base_cgf$rsim(
        n = n,
        vector_length = vector_length,
        parameter_vector = parameter_vector,
        tvec = t_total,
        flatten = FALSE,
        ...
      )
    }
  }

  # call_history label
  hist <- paste(base_cgf$call_history, collapse = " -> ")
  op_name <- c(hist, "ExponentialTiltCGF")

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
    K2operatorAK2AT = K2operatorAK2AT,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    ineq_constraint = ineq_constraint,
    analytic_tvec_hat = analytic_tvec_hat,
    rsim = rsim,
    op_name = op_name
  )

  do.call(createCGF, c(cgf_args, list(...)))
}

#' @title Exponential Tilting of a CGF
#'
#' @description
#' Returns a new `CGF` corresponding to the exponential tilt of a base CGF.
#'
#' If the base random vector has CGF \eqn{K(t;\theta)}, then the exponentially tilted CGF is
#' \deqn{K_{\text{tilt}}(t;\theta) = K(t + h(\theta);\theta) - K(h(\theta);\theta).}
#'
# This is not a randomly-stopped sum (RSS). It does not introduce a random count;
# it changes the measure by exponential reweighting.
#'
#' @param cgf A `CGF` object.
#' @param tilt A numeric vector, an `adaptor`, or a function mapping `theta -> h`.
#'   - If numeric, it is treated as a fixed tilt vector.
#'   - If an adaptor/function, it must return the tilt vector \eqn{h(\theta)}.
#'   The returned `h` must have length 1, or divide `length(tvec)` at evaluation time.
#' @param iidReps Either `"any"` or a positive integer.
#' @param block_size Optional block size for iid replication.
#' @param ... Passed to `createCGF` (advanced overrides).
#'
#' @return A `CGF` object (tilted).
#'
#' @examples
#' ## ------------------------------------------------------------
#' ## Example 1: Poisson — tilt gives another Poisson (exact)
#' ## ------------------------------------------------------------
#' lam <- 2
#' h   <- 0.4
#' cgT <- ExponentialTiltCGF(PoissonCGF, tilt = h, iidReps = 1)
#' tt <- c(-0.2, 0.1, 0.05)
#' stopifnot(all.equal(cgT$K(tt, lam), PoissonCGF$K(tt, lam * exp(h)), tol = 1e-12))
#'
#' ## ------------------------------------------------------------
#' ## Example 2: Gamma — tilt shifts the rate (exact)
#' ##   Gamma(a,b): K(tt) = -a log(1 - tt/b), domain tt < b
#' ##   Tilt by h: Gamma(a, b-h), domain tt < (b-h)
#' ## ------------------------------------------------------------
#' a <- 5; b <- 3; h <- 0.5
#' cgTg <- ExponentialTiltCGF(GammaCGF, tilt = h, iidReps = 1)
#' stopifnot(all.equal(cgTg$K(0.2, c(a,b)), GammaCGF$K(0.2, c(a, b - h)), tol = 1e-12))
#' # Check constraint: requires both tt+h < b and h < b
#' print(cgTg$ineq_constraint(0.2, c(a,b)))  # should be <= 0
#'
#' ## ------------------------------------------------------------
#' ## Example 3: No closed-form family - tilt a sum CGF
#' ##   Z = Pois(lambda) + Binom(n,p)  (convolution, not a named family)
#' ##   But we can still compute the exact likelihood numerically.
#' ## ------------------------------------------------------------
#' set.seed(1)
#' n_fix <- 10
#' cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
#' cg_bin  <- BinomialModelCGF(n = adaptor(fixed_param = n_fix),
#'                             prob = adaptor(indices = 2), iidReps = 1)
#' cg_sum  <- sumOfIndependentCGF(list(cg_pois, cg_bin), iidReps = 1)
#'
#' h <- 0.3
#' cg_sum_tilt <- ExponentialTiltCGF(cg_sum, tilt = h, iidReps = "any", block_size = 1)
#'
#' # Simulate from the exact tilted model:
#' theta_base <- c(lambda = 2.5, p = 0.35)
#' lambda_t <- theta_base[1] * exp(h)
#' p_t      <- plogis(qlogis(theta_base[2]) + h)
#' B <- 50
#' y <- rpois(B, lambda_t) + rbinom(B, n_fix, p_t)
#'
#' # Exact likelihood via convolution (finite sum):
#' dZ <- function(z, lambda, p) {
#'   k <- 0:n_fix
#'   sum(dbinom(k, n_fix, p) * dpois(z - k, lambda))
#' }
#' nll_exact <- function(theta) {
#'   lam <- theta[1]; p <- theta[2]
#'   lam_t <- lam * exp(h)
#'   p_t   <- plogis(qlogis(p) + h)
#'   -sum(log(vapply(y, dZ, numeric(1), lambda = lam_t, p = p_t)))
#' }
#'
#' # SPA MLE using the tilted CGF
#' cgB <- iidReplicatesCGF(cg_sum_tilt, iidReps = B, block_size = 1)
#' fit_spa <- find.saddlepoint.MLE(
#'   observed.data  = y,
#'   cgf            = cgB,
#'   starting.theta = c(1, 0.5),
#'   lb.theta       = c(1e-8, 1e-8),
#'   ub.theta       = c(Inf, 1-1e-8),
#'   method         = "two_step"
#' )
#'
#' # Exact MLE (numeric) for comparison:
#' fit_exact <- nlminb(
#'   start = c(1, 0.5),
#'   objective = nll_exact,
#'   lower = c(1e-8, 1e-8),
#'   upper = c(Inf, 1 - 1e-8)
#' )
#'
#' cat("True theta_base:", theta_base, "\n")
#' cat("SPA  theta_hat :", fit_spa$MLEs.theta, "\n")
#' cat("Exact theta_hat:", fit_exact$par, "\n")
#'
#' ## ------------------------------------------------------------
#' ## Example 4: theta-dependent tilt + SPA MLE vs exact MLE
#' ##   Z = Pois(lambda) + Binom(n,p), then exponential-tilt by h(theta)
#' ##   where h(theta) = log(lambda).
#' ##
#' ##   Under tilt by h, the components tilt as:
#' ##     Pois(lambda)    -> Pois(lambda * exp(h))
#' ##     Binom(n,p)      -> Binom(n, plogis(qlogis(p) + h))
#' ##
#' ##   Here h depends on theta, so the tilted likelihood is still a
#' ##   well-defined model in theta, but it is NOT “just a fixed tilted
#' ##   family” anymore.
#' ## ------------------------------------------------------------
#' \donttest{
#' set.seed(2)
#' B     <- 60
#' n_fix <- 10
#'
#' ## Base CGF for one observation Z = X + Y
#' cg_pois <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
#' cg_bin  <- BinomialModelCGF(
#'   n    = adaptor(fixed_param = n_fix),
#'   prob = adaptor(indices = 2),
#'   iidReps = 1
#' )
#' cg_sum <- sumOfIndependentCGF(list(cg_pois, cg_bin), iidReps = 1)
#'
#' ## Theta-dependent tilt:
#' ## h(theta) = log(lambda).  (Requires lambda>0; we enforce via lb.theta.)
#' tilt_theta <- function(theta) log(theta[1])
#'
#' ## tilted CGF for B i.i.d. observations (scalar blocks)
#' cg_tilt_B <- ExponentialTiltCGF(
#'   cgf        = cg_sum,
#'   tilt       = tilt_theta,
#'   iidReps    = B,
#'   block_size = 1
#' )
#'
#' ## Simulate data from the exact tilted model at theta_true:
#' theta_true <- c(lambda = 2.5, p = 0.35)
#' h_true     <- tilt_theta(theta_true)
#'
#' lambda_t <- theta_true[1] * exp(h_true)                 # = lambda^2
#' p_t      <- plogis(qlogis(theta_true[2]) + h_true)
#'
#' y <- rpois(B, lambda_t) + rbinom(B, n_fix, p_t)
#'
#' ## SPA MLE under the tilted CGF
#' fit_spa <- find.saddlepoint.MLE(
#'   observed.data  = y,
#'   cgf            = cg_tilt_B,
#'   starting.theta = c(1.5, 0.5),
#'   lb.theta       = c(1e-8, 1e-8),
#'   ub.theta       = c(Inf, 1 - 1e-8),
#'   method         = "two_step"
#' )
#'
#' ## finite convolution for comparison
#' dZ <- function(z, lambda, p) {
#'   k <- 0:n_fix
#'   sum(dbinom(k, n_fix, p) * dpois(z - k, lambda))
#' }
#'
#' nll_exact <- function(theta) {
#'   h       <- tilt_theta(theta)
#'   lam_t   <- theta[1] * exp(h)
#'   p_tilt  <- plogis(qlogis(theta[2]) + h)
#'   likvec  <- vapply(y, dZ, numeric(1), lambda = lam_t, p = p_tilt)
#'   -sum(log(pmax(likvec, .Machine$double.xmin)))
#' }
#'
#' fit_exact <- nlminb(
#'   start     = c(1.5, 0.5),
#'   objective = nll_exact,
#'   lower     = c(1e-8, 1e-8),
#'   upper     = c(Inf, 1 - 1e-8)
#' )
#'
#' cat("theta_true        =", paste(round(theta_true, 4), collapse=" "), "\n")
#' cat("theta_hat_spa     =", paste(round(fit_spa$MLEs.theta, 4), collapse=" "), "\n")
#' cat("theta_hat_exact   =", paste(round(fit_exact$par, 4), collapse=" "), "\n")
#' }
#'
#' @export
ExponentialTiltCGF <- function(cgf,
                           tilt,
                           iidReps = "any",
                           block_size = NULL,
                           ...) {
  stopifnot(inherits(cgf, "CGF"))
  .check_iidReps(iidReps)

  # Allow numeric tilt directly (fixed tilt)
  if (is.numeric(tilt) && !is.function(tilt)) {
    tilt <- adaptor(fixed_param = as.numeric(tilt))
  }

  tilt_fn <- validate_function_or_adaptor(tilt)

  base <- .exponentialTiltCGF_internal(base_cgf = cgf, tilt_fn = tilt_fn, ...)
  iidReplicatesCGF(cgf = base, iidReps = iidReps, block_size = block_size)
}
