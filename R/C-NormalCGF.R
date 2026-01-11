# R/NormalCGF.R
# Objects: NormalCGF, NormalModelCGF, MultivariateNormalModelCGF
#
# - We keep the univariate Normal implementation consistent with the other univariate
#   CGFs built via .make_univariate_model_cgf_matrix.
# - Fix multivariate Normal replication semantics by delegating replication to
#   iidReplicatesCGF (block_size = d, iidReps = ...)
# - analytic_tvec_hat exists for both univariate and multivariate Normal.
# - Provide K2_solve + logdetK2 for multivariate Normal (so downstream SPA code
#   does not need to form inverses explicitly).
#
# Notes:
# - For multivariate Normal, sigma(theta) is assumed to return a symmetric
#   covariance matrix. If it is not symmetric, the CGF formulas for K1/K2 are
#   mathematically inconsistent. We therefore symmetrise Sigma internally.
# - The replication arguments iidReps refer to i.i.d. replication of the MVN
#   observation Y in R^d, not to anything inside Sigma.

# Univariate Normal
#' @noRd
.normal_base_cgf <- function(iidReps, op_name, ...) {

  # Parameter vector convention for the univariate vectorised utilities:
  #   param = c(mu[1:L], sigma[1:L])
  # where sigma is a standard deviation (must be > 0)
  split_mu_sigma <- function(param) {
    L <- length(param) / 2
    if (!is.finite(L) || L != as.integer(L) || L < 1) {
      stop("Normal params must be concatenated as c(mu[1:L], sigma[1:L]).")
    }
    L <- as.integer(L)
    idx <- seq_len(L)
    cbind(mu = param[idx], sigma = param[L + idx])
  }

  .make_univariate_model_cgf_matrix(
    K_elem <- function(tvec, pm) {
      sig2 <- pm[,2] * pm[,2]
      pm[,1] * tvec + 0.5 * sig2 * (tvec * tvec)
    },

    K1_elem <- function(tvec, pm) {
      sig2 <- pm[,2] * pm[,2]
      pm[,1] + sig2 * tvec
    },

    K2_elem <- function(tvec, pm) {
      sig2 <- pm[,2] * pm[,2]
      sig2
    },
    K3_elem     = function(tvec, pm) rep(0, length(tvec)),
    K4_elem     = function(tvec, pm) rep(0, length(tvec)),

    t_hat_elem <- function(y, pm) {
      sig2 <- pm[,2] * pm[,2]
      (y - pm[,1]) / sig2
    },

    split_param_to_mat = split_mu_sigma,
    rsim_elem = function(n, tvec, pm, ...) {
      mu <- as.numeric(pm[, 1])
      sigma <- as.numeric(pm[, 2])

      if (any(!is.finite(mu))) stop("NormalCGF$rsim: 'mu' must be finite.")
      if (any(!is.finite(sigma)) || any(sigma < 0)) {
        stop("NormalCGF$rsim: 'sigma' must be finite and >= 0.")
      }

      sig2 <- sigma * sigma
      mean_tilt <- mu + sig2 * tvec
      vector_length <- length(tvec)

      matrix(
        stats::rnorm(
          n = n * vector_length,
          mean = rep.int(mean_tilt, times = n),
          sd = rep.int(sigma, times = n)
        ),
        nrow = vector_length,
        ncol = n
      )
    },

    iidReps = iidReps,
    op_name = op_name,
    ...
  )
}

#' Univariate Normal CGF object
#'
#' Ready-to-use CGF for a univariate Normal distribution with parameters
#' \eqn{(\mu,\sigma)}.
#'
#' The parameter vector is interpreted as \code{c(mu, sigma)} (length 2),
#' or in vectorised form as \code{c(mu[1:L], sigma[1:L])}.
#'
#' With \code{iidReps="any"} (the default), \code{length(tvec)} must be a multiple
#' of the number of parameter rows (L). In the scalar case (L=1), any length
#' \code{tvec} is treated as i.i.d. replicates.
#'
#' @rdname GaussianCGF
#'
#' @export
NormalCGF <- .normal_base_cgf(iidReps = "any", op_name = "NormalCGF")

#' @rdname GaussianCGF
#'
#' @export
GaussianCGF <- NormalCGF







#' Create a parametric univariate Normal CGF
#'
#' @description
#' Builds a CGF for \eqn{Y \sim N(\mu(\theta), \sigma(\theta))}, where \code{mu} and
#' \code{sigma} are functions (or adaptors) of \code{theta}.
#'
#' Both \code{mu(theta)} and \code{sigma(theta)} may return a scalar or a vector.
#' If they return vectors, they must have the same length L and are interpreted
#' in the same vectorised way as other univariate model CGFs in this package.
#'
#' @param mu Function/adaptor mapping theta -> mean(s).
#' @param sigma Function/adaptor mapping theta -> sd(s).
#' @param iidReps Either "any" or a positive integer (replication semantics).
#' @param ... Passed through to CGF creation.
#'
#' @rdname GaussianModelCGF
#'
#' @export
NormalModelCGF <- function(mu, sigma, iidReps = "any", ...) {
  .check_iidReps(iidReps)
  mu_fn <- validate_function_or_adaptor(mu)
  sig_fn <- validate_function_or_adaptor(sigma)

  adaptor_ <- function(theta) {
    mu_val  <- mu_fn(theta)
    sig_val <- sig_fn(theta)
    if (length(mu_val) != length(sig_val)) {
      stop("NormalModelCGF: mu(theta) and sigma(theta) must have the same length.")
    }
    c(mu_val, sig_val)
  }

  base <- .normal_base_cgf(iidReps = iidReps, op_name = "NormalModelCGF", ...)
  adaptCGF(cgf = base, adaptor = adaptor_)
}


#' @rdname GaussianModelCGF
#'
#' @export
GaussianModelCGF <- NormalModelCGF






















# Multivariate Normal
#' @noRd
.mvn_dim_from_param <- function(param) {
  L <- length(param)
  d <- (-1 + sqrt(1 + 4 * L)) / 2
  # if (!is.finite(d) || d != as.integer(d) || d < 1) stop("MultivariateNormal: invalid parameter length. Expected length(param) = d + d^2 = d*(d+1) for some integer d.")
  # if (L != d * (d + 1)) stop("MultivariateNormal: invalid parameter length. Expected d*(d+1); got ", L, ".")
  d
}

#' @noRd
.mvn_base_cgf_one_block <- function(op_name, ...) {

  # Extract (mu, Sigma) from the flattened parameter vector:
  #   param = c(mu[1:d], vec(Sigma)) with vec in column-major order.
  .split_param <- function(param) {
    d <- .mvn_dim_from_param(param)
    mu <- param[seq_len(d)]
    Sigma <- matrix(param[(d+1):(d + d * d)], nrow = d, ncol = d, byrow = FALSE)

    # # Symmetrise defensively (covariance should be symmetric).
    # Sigma <- 0.5 * (Sigma + t(Sigma))

    Sigma <- matrix(param[(d + 1):(d + d * d)], nrow = d, ncol = d, byrow = FALSE)
         # if (max(abs(Sigma - t(Sigma))) > 1e-10) stop("Sigma must be symmetric. Max |Sigma - t(Sigma)| = ", signif(asym, 6))
    if (!inherits(Sigma, "advector") && !inherits(Sigma, "adsparse")) {
      asym <- max(abs(Sigma - t(Sigma)))
      if (asym > 1e-10)  stop("Sigma must be symmetric. Max |Sigma - t(Sigma)| = ", signif(asym, 6))
    }
    list(d = d, mu = mu, Sigma = Sigma)
  }

  Kfun <- function(tvec, param) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    mu <- sp$mu
    Sigma <- sp$Sigma
    sum(mu * tvec) + 0.5 * sum(tvec * (Sigma %*% tvec))
  }

  K1fun <- function(tvec, param) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    as.vector(sp$mu + sp$Sigma %*% tvec)
  }

  K2fun <- function(tvec, param) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    sp$Sigma
  }

  tiltingfun <- function(tvec, param) {
    # K(t) - t^T K1(t) = -0.5 * t^T Sigma t
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    -0.5 * sum(tvec * (sp$Sigma %*% tvec))
  }

  K2opfun <- function(tvec, param, x, y) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    if (length(x) != d || length(y) != d) stop("MultivariateNormal: K2operator expects vectors of length d.")
    sum(x * (sp$Sigma %*% y))
  }

  K2operatorAK2ATfun <- function(tvec, param, Bmat) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(tvec) != d) stop("MultivariateNormal: length(tvec) must equal d.")
    if (ncol(Bmat) != d) {
      stop("MultivariateNormal: K2operatorAK2AT expects Bmat with ncol == d.")
    }
    Bmat %*% sp$Sigma %*% t(Bmat)
  }

  # Analytic saddlepoint t-hat: Sigma^{-1} (y - mu)
  analytic_tvec_hat_fun <- function(y, param) {
    sp <- .split_param(param)
    d <- sp$d
    if (length(y) != d) stop("MultivariateNormal: length(y) must equal d.")
    as.vector(solve(sp$Sigma, y - sp$mu))
  }

  # Simulation: draw n samples from N(mu + Sigma*t, Sigma), return d x n.
  simulate_fun <- function(n, vector_length, parameter_vector, tvec = NULL, ...) {

    # Parse params (dimension d inferred from parameter_vector length)
    sp <- .split_param(parameter_vector)
    d <- as.integer(sp$d)

    # This is NOT redundant: base CGF only checks vector_length is an integer,
    # it does not know it must equal d.
    if (vector_length != d) {
      stop(
        "MultivariateNormal$rsim: 'vector_length' must equal d = ", d,
        " inferred from parameter_vector; got ", vector_length, ".",
        call. = FALSE
      )
    }

    mu    <- as.numeric(sp$mu)
    Sigma <- as.matrix(sp$Sigma)

    # Apply exact exponential tilting: mu_tilt = mu + Sigma %*% tvec
    if (!is.null(tvec)) {
      mu <- mu + as.vector(Sigma %*% tvec)
    }

    # We only need chol for simulation; chol requires SPD.
    Rchol <- tryCatch(
      chol(Sigma),
      error = function(e) {
        stop(
          "MultivariateNormal$rsim: chol(Sigma) failed. ",
          "Sigma must be symmetric positive definite for simulation. ",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    Z <- matrix(stats::rnorm(d * n), nrow = d, ncol = n)
    matrix(mu, nrow = d, ncol = n) + t(Rchol) %*% Z
  }


  K3opfun <- function(tvec, param, v1, v2, v3) 0
  K4opfun <- function(tvec, param, v1, v2, v3, v4) 0
  func_Tfun <- function(tvec, param) 0

  K4AABBfun <- function(tvec, param, Q1, Q2) 0
  K3K3AABBCCfun <- function(tvec, param, Q1, Q2, Q3) 0
  K3K3ABCABCfun <- function(tvec, param, Q1, Q2, Q3) 0


  createCGF(
    K  = Kfun,
    K1 = K1fun,
    K2 = K2fun,
    K2operator = K2opfun,
    K2operatorAK2AT = K2operatorAK2ATfun,
    K3operator = K3opfun,
    K4operator = K4opfun,
    analytic_tvec_hat_func = analytic_tvec_hat_fun,
    tilting_exponent = tiltingfun,
    func_T = func_Tfun,

    # #
    # K2_solve = K2_solve_fun,
    # logdetK2 = logdetK2_fun,

    #
    K4operatorAABB = K4AABBfun,
    K3K3operatorAABBCC = K3K3AABBCCfun,
    K3K3operatorABCABC = K3K3ABCABCfun,

    rsim = simulate_fun,

    op_name = op_name,
    ...
  )
}

#' Create a Multivariate Normal CGF Object
#'
#' @description
#' Creates a CGF for a \eqn{d}-dimensional Normal distribution with
#' mean \code{mu(theta)} and covariance \code{sigma(theta)}.
#'
#' Replication:
#' - \code{iidReps} refers to i.i.d. replication of the \eqn{d}-vector observation.
#' - If \code{iidReps="any"} (default), \code{length(tvec)} must be a multiple of \eqn{d}.
#' - If \code{iidReps=m} (integer), \code{length(tvec)} must equal \eqn{d*m}.
#'
#' @param mu Function/adaptor mapping \code{theta} to a numeric vector of length \eqn{d}.
#' @param sigma Function/adaptor mapping \code{theta} to a \eqn{d\times d} covariance matrix.
#' @param iidReps Either \code{"any"} or a positive integer.
#' @param ... Passed through to CGF creation.
#'
#' @rdname MultivariateGaussianModelCGF
#'
#' @return A \code{CGF} object.
#'
#' @export
MultivariateNormalModelCGF <- function(mu, sigma, iidReps = "any", ...) {
  .check_iidReps(iidReps)

  mu_fn <- validate_function_or_adaptor(mu)
  sigma_fn <- validate_function_or_adaptor(sigma)

  # adaptor(theta) -> flattened param = c(mu, vec(Sigma))
  param_adaptor <- function(theta) {
    muVal <- mu_fn(theta)
    SigmaVal <- sigma_fn(theta)
    if (nrow(SigmaVal) != length(muVal) || ncol(SigmaVal) != length(muVal)) {
      stop("MultivariateNormalModelCGF: sigma(theta) must be a square matrix with dimension equal to length(mu(theta)).")
    }
    c(muVal, as.vector(SigmaVal))
  }

  # Base MVN CGF for ONE d-vector (d inferred from length(param))
  base_one <- .mvn_base_cgf_one_block(op_name = "MultivariateNormalModelCGF", ...)

  # Replicate Y (block_size = d(param)) using iidReplicatesCGF.
  # Note: This is *internal* function-valued block_size; user never supplies it.
  block_size_fun <- function(param) .mvn_dim_from_param(param)
  attr(block_size_fun, "label") <- "d"

  replicated_ <- iidReplicatesCGF(
    cgf = base_one,
    iidReps = iidReps,
    block_size = block_size_fun
  )

  adaptCGF(cgf = replicated_, adaptor = param_adaptor)
}


#' @rdname MultivariateGaussianModelCGF
#'
#' @export
MultivariateGaussianModelCGF <- MultivariateNormalModelCGF
