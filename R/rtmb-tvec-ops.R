# R/rtmb-tvec-ops.R
# Utilities to obtain t-hat(theta) under RTMB, with AD captured.


.rtmb_spa_env <- local({
  e <- new.env(parent = emptyenv())
  e$H_user   <- list()  # optional cache slots if you decide to key / reuse
  e$H_solver <- list()
  e$y        <- NULL    # streamed into tapes via DataEval
  e
})
rtmb_spa_cache_clear <- function() {
  rm(list = ls(.rtmb_spa_env), envir = .rtmb_spa_env)
  .rtmb_spa_env$H_user   <- list()
  .rtmb_spa_env$H_solver <- list()
  .rtmb_spa_env$y        <- NULL
  invisible(TRUE)
}





#' Solve (K2)^{-1} w at (t, theta)
#' @noRd
create_tvec_hat_K2_solve_fn <- function(cgf) {
  force(cgf)
  function(tvec, theta, w) {
    # RTMB knows solve() for AD types; see AD matrix methods. 
    # Use as-is: this is AD-safe when tvec/theta/w are AD.  (RTMB manual: ADmatrix) 
    solve(cgf$K2(tvec, theta), w)
  }
}




# #' @noRd
# build_tvec_hat_newton_tape <- function(cgf, y, theta_init, t_init = NULL, ...)
# {
#   stopifnot(inherits(cgf, "CGF"), is.numeric(theta_init), is.numeric(y))
#   m <- length(y)
#   if (is.null(t_init)) t_init <- rep(0, m)
# 
#   # Stream y from the environment (can be updated between evaluations)
#   .rtmb_spa_env$y <- as.numeric(y)
# 
#   # Objective in (t, theta); grad_wrt_t = K1 - y
#   obj_fun <- function(z) {
#     u     <- z[seq_len(m)]
#     theta <- z[-seq_len(m)]
#     # Anchor DataEval to the tape by passing an AD argument (u)
#     yAD <- RTMB::DataEval(function(dummy) .rtmb_spa_env$y, u)
#     # yAD <- RTMB::DataEval(function(i) .rtmb_spa_env$y[i], seq_len(m))
#     cgf$K(u, theta) - sum(u * yAD)
#   }
# 
#   # Tape in (t, theta)
#   F_tape <- RTMB::MakeTape(obj_fun, x = c(t_init, theta_init))
# 
# 
#   Uhat_theta <- F_tape$newton(random = seq_len(m), ...)
# 
#   # Wrap to a pure theta → t̂(theta) tape
#   G <- RTMB::MakeTape(function(theta) Uhat_theta(theta), theta_init)
#   attr(G, "t_init") <- t_init
#   G
# }




#' @noRd
.make_strictly_feasible_t0 <- function(cgf, t_init, theta, interior_margin = 1e-10) {
  g <- cgf$ineq_constraint(t_init, theta)
  if (!length(g)) return(t_init)                      # no constraints
  if (all(is.finite(g)) && max(g) < -interior_margin) return(t_init)
  
  base <- rep(0, length(t_init))  # 0 is interior for Gamma/Exp and their linear maps
  step <- 1.0
  for (k in 1:60) {
    t_try <- (1 - step) * base + step * t_init
    g_try <- cgf$ineq_constraint(t_try, theta)
    if (all(is.finite(g_try)) && max(g_try) < -interior_margin) return(t_try)
    step <- step / 2
  }
  stop("Could not find a strictly feasible starting t (log-barrier needs g(t,theta)<0).")
}

#' @noRd
auto_mu_from_grad_balance <- function(cgf, y, theta, t0, eta = 0.05) {
  gvals0 <- cgf$ineq_constraint(t0, theta)
  if (!length(gvals0)) return(0)                           # no constraints
  if (!all(is.finite(gvals0)) || any(gvals0 >= 0))
    stop("t0 must be strictly feasible for the log barrier (g(t0,theta) < 0).")
  
  # Main gradient size
  g_main <- cgf$K1(t0, theta) - y
  n_main <- sqrt(sum(g_main^2))
  if (!is.finite(n_main) || n_main == 0) n_main <- 1
  
  # Barrier gradient size at start: -sum_i (1/g_i) * grad g_i
  m <- length(t0)
  ineq_tape <- RTMB::MakeTape(
    function(a) cgf$ineq_constraint(a[seq_len(m)], a[-seq_len(m)]),
    x = c(t0, theta)
  )
  J <- ineq_tape$jacobian(c(t0, theta))[, seq_len(m), drop = FALSE]  # grad_g/grat_t as rows
  grad_phi <- -colSums(J / as.numeric(gvals0))  # broadcast by rows
  n_phi <- sqrt(sum(grad_phi^2))
  if (!is.finite(n_phi) || n_phi == 0) return(0)
  
  mu <- eta * n_main / n_phi
  max(min(mu, 1e+2), 1e-8)   # conservative clamps
}

#' @noRd
build_tvec_hat_newton_tape <- function(
    cgf, y, theta_init, t_init = NULL,
    barrier = c("log","soft","none"),
    mu = "auto", eta = 0.05, interior_margin = 1e-10, ...
) {
  barrier <- match.arg(barrier)
  m <- length(y)
  if (is.null(t_init)) t_init <- rep(0, m)
  
  # stream data (y, μ, τ) into the tape
  .rtmb_spa_env$y <- as.numeric(y)
  
  has_ineq <- function(t_vec, th) length(cgf$ineq_constraint(t_vec, th)) > 0L
  
  # If no constraints exist, override to "none"
  if (!has_ineq(t_init, theta_init) && barrier != "none") barrier <- "none"
  
  # Feasible start for log barrier
  t0 <- if (barrier == "log") .make_strictly_feasible_t0(cgf, t_init, theta_init, interior_margin) else t_init
  
  # Choose mu and (if needed) tau outside the tape; stream with DataEval
  if (barrier == "log") {
    mu_val <- if (identical(mu, "auto")) auto_mu_from_grad_balance(cgf, y, theta_init, t0, eta) else as.numeric(mu)
    .rtmb_spa_env$mu <- mu_val
  } else if (barrier == "soft") {
    g0 <- cgf$ineq_constraint(t0, theta_init)
    if (!length(g0)) { .rtmb_spa_env$tau <- 1; .rtmb_spa_env$mu <- 0 } else {
      s0 <- max(-g0)                         # initial slack to boundary
      tau <- max(1e-12, 0.05 * s0)
      .rtmb_spa_env$tau <- tau
      .rtmb_spa_env$mu  <- 1.0
    }
  } else {
    .rtmb_spa_env$mu <- 0
  }
  
  # Smooth barrier term added to the (t,theta) objective; y,mu,tau flow via DataEval
  barrier_term <- function(t_vec, theta) {
    # NOTE: no comparisons on AD types anywhere in this function.
    # 'barrier' is a plain R string (non-AD), set outside the tape.
    g <- cgf$ineq_constraint(t_vec, theta)
    if (!length(g)) return(0)   # length() is metadata, not AD
    
    muAD <- RTMB::DataEval(function(dummy) .rtmb_spa_env$mu, t_vec)
    
    if (barrier == "log") {
      # log barrier: requires strict feasibility; we enforce via t0 outside the tape
      # no check of 'muAD' here; just multiply
      -muAD * sum(log(-g))
    } else if (barrier == "soft") {
      # softplus barrier: fully smooth and defined everywhere
      tauAD <- RTMB::DataEval(function(dummy) .rtmb_spa_env$tau, t_vec)
      muAD * sum(log1p(exp(g / tauAD)))
    } else {
      # "none"
      0
    }
  }
  
  
  # Objective in (t, theta); grad_wrt_t is K1 - y + barrier_grad
  obj_fun <- function(z) {
    u     <- z[seq_len(m)]
    theta <- z[-seq_len(m)]
    yAD   <- RTMB::DataEval(function(dummy) .rtmb_spa_env$y, u)
    cgf$K(u, theta) - sum(u * yAD) + barrier_term(u, theta)
  }
  
  # Tape over (t,theta)
  F_tape <- RTMB::MakeTape(obj_fun, x = c(t0, theta_init))
  # Optimize out t with Newton, producing a θ -> t̂(θ) function we can tape again
  Uhat_theta <- F_tape$newton(random = seq_len(m), ...)
  
  # Wrap to a pure theta → t̂(theta) tape (so outer optimizers see a clean mapping)
  G <- RTMB::MakeTape(function(theta) Uhat_theta(theta), theta_init)
  attr(G, "t_init")  <- t0
  attr(G, "barrier") <- barrier
  attr(G, "mu")      <- get0("mu", envir = .rtmb_spa_env, inherits = FALSE)
  G
}













# ===== Atomics for the other two modes =======================================

# -----------------------------------------------------------------------------
# Fixed user-supplied tvec0: Hessian-safe ADjoint
#   forward  : f(theta) = tvec0  (user value; one-time tape recommended)
#   reverse  : -(∂_θ K1(y,θ))^T * w,  where  w = K2(y,θ)^{-1} * dy  and y is tape-time t
# Implementation detail:
#   Pre-tape Htheta(theta, t, w) = ∂_θ [ w' * K1(t, theta) ] and embed atomically.
#   This keeps dependence on (theta, t, w) so higher orders are correct.
# -----------------------------------------------------------------------------
#' @noRd
make_user_tvec_atomic <- function(tvec0, cgf, theta_init) {
  stopifnot(is.numeric(tvec0), is.numeric(theta_init))
  p <- length(theta_init)
  m <- length(tvec0)
  
  # Pre-tape H_scalar(theta, t, w) = w' * K1(t, theta)
  H_scalar <- RTMB::MakeTape(
    function(x) {
      theta <- x[seq_len(p)]
      tvec  <- x[p + seq_len(m)]
      w     <- x[p + m + seq_len(m)]
      sum(w * cgf$K1(tvec, theta))
    },
    x = c(theta_init, tvec0, rep(0, m))
  )
  Hgrad  <- H_scalar$jacfun()  # gradient wrt (theta, t, w); returns length p+m+m
  Htheta <- RTMB::MakeTape(
    function(x) {
      g <- Hgrad(x)
      g[seq_len(p)]           # extract the theta-block; still a *tape* in (theta, t, w)
    },
    x = c(theta_init, tvec0, rep(0, m))
  )
  Htheta_atm <- Htheta$atomic()  # embeddable, higher-order safe
  
  # Forward: the user-supplied t-value (one-time tape semantics)
  f  <- function(theta) tvec0
  
  # Reverse: evaluate at the *tape-time* y (not frozen tvec0)
  df <- function(theta, y, dy) {
    thetaAD <- RTMB::AD(theta)
    yAD     <- RTMB::AD(y)
    dyAD    <- RTMB::AD(dy)
    
    # AD-aware linear solve: w = K2(y,theta)^{-1} * dy
    wAD <- solve(cgf$K2(yAD, thetaAD), dyAD)
    
    # -(∂_θ K1(y,θ))^T w  via the taped gradient (depends on theta, t=y, w)
    -as.vector( Htheta_atm(c(thetaAD, yAD, wAD)) )
  }
  
  RTMB::ADjoint(f, df, name = "user_tvec_atomic")
}









#' @noRd
make_solver_tvec_atomic <- function(cgf, y, solver_fun, theta_init, t_init = NULL) {
  m <- length(y)
  if (is.null(t_init)) t_init <- rep(0, m)
  p <- length(theta_init)
  
  # Pre-tape H_scalar(theta, t, w) = w' K1(t, theta)
  H_scalar <- RTMB::MakeTape(
    function(x) {
      theta <- x[seq_len(p)]
      tvec  <- x[p + seq_len(m)]
      w     <- x[p + m + seq_len(m)]
      sum(w * cgf$K1(tvec, theta))
    },
    x = c(theta_init, t_init, rep(0, m))
  )
  Hgrad      <- H_scalar$jacfun()  # gradient wrt (theta, t, w)
  Htheta     <- RTMB::MakeTape(function(x) Hgrad(x)[seq_len(p)],
                               x = c(theta_init, t_init, rep(0, m)))
  Htheta_atm <- Htheta$atomic()
  
  # Forward: call the numeric solver each time
  f  <- function(theta) {
    solver_fun(theta, y = y, cgf = cgf, starting.tvec = t_init)
  }
  
  # Reverse: implicit-function formula using AD primitives only
  df <- function(theta, tvec, dtvec) {
    thetaAD <- RTMB::AD(theta)
    tvecAD  <- RTMB::AD(tvec)
    dtAD    <- RTMB::AD(dtvec)
    wAD     <- solve(cgf$K2(tvecAD, thetaAD), dtAD)
    -as.vector(Htheta_atm(c(thetaAD, tvecAD, wAD)))
  }
  
  RTMB::ADjoint(f, df, name = "solver_tvec_atomic")
}






