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

#' Build a Tape mapping theta -> t_hat(theta) via Newton on
#'   F_tape(t,theta) = K(t;theta) - sum(t * y).
#' Returns an RTMB Tape 'G' with G(theta) == t_hat(theta).
#' @noRd

# build_tvec_hat_newton_tape <- function(cgf, y, theta_init, t_init = NULL, ...) {
#   stopifnot(inherits(cgf, "CGF"), is.numeric(theta_init), is.numeric(y))
#   m <- length(y)
#   if (is.null(t_init)) t_init <- rep(0, m)
#   
#   # Scalar objective whose gradient wrt t is K1 - y. Convex in t for CGFs.
#   obj_fun <- function(z) {
#     u <- z[seq_len(m)]
#     theta <- z[-seq_len(m)]
#     cgf$K(u, theta) - sum(u * y)
#   }
#   
#  
#   
#   # Tape for the joint objective in (t, theta), initialized at (t_init, theta_init)
#   F_tape <- RTMB::MakeTape(obj_fun, x = c(t_init, theta_init))
#   
#   # NOTE: In RTMB >= 1.7 the argument name is 'random' (not 'indices')
#   # Returns a tape Uhat_theta : theta -> t_hat(theta)
#   Uhat_theta <- F_tape$newton(random = seq_len(m), ...)
#   
#   # Wrap to ensure a clean theta -> t_hat(theta) map (no concatenation with t_init)
#   G <- RTMB::MakeTape(function(theta) Uhat_theta(theta), theta_init)
#   
#   attr(G, "t_init") <- t_init
#   G
# }

build_tvec_hat_newton_tape <- function(cgf, y, theta_init, t_init = NULL, ...) {
  stopifnot(inherits(cgf, "CGF"), is.numeric(theta_init), is.numeric(y))
  m <- length(y)
  if (is.null(t_init)) t_init <- rep(0, m)
  
  # Stream y from the environment (can be updated between evaluations)
  .rtmb_spa_env$y <- as.numeric(y)
  
  # Objective in (t, theta); grad_wrt_t = K1 - y
  obj_fun <- function(z) {
    u     <- z[seq_len(m)]
    theta <- z[-seq_len(m)]
    # Anchor DataEval to the tape by passing an AD argument (u)
    yAD <- RTMB::DataEval(function(dummy) .rtmb_spa_env$y, u)
    # yAD <- RTMB::DataEval(function(i) .rtmb_spa_env$y[i], seq_len(m))
    cgf$K(u, theta) - sum(u * yAD)
  }
  
  # Tape in (t, theta)
  F_tape <- RTMB::MakeTape(obj_fun, x = c(t_init, theta_init))
  
  # RTMB ≥ 1.7 uses argument 'random' (higher‑order‑safe Newton)
  Uhat_theta <- F_tape$newton(random = seq_len(m), ...)
  
  # Wrap to a pure theta → t̂(theta) tape
  G <- RTMB::MakeTape(function(theta) Uhat_theta(theta), theta_init)
  attr(G, "t_init") <- t_init
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






