# R/A-CGF_fromMultinomialFamily.R
# Objects: MultinomialFamilyCGF, createMultinomialFamilyCGF



# Strictly speaking, all we need for the CGF is of the form
# K(t) = N log ( sum_i p_i e^{t_i} )
# If we parameterize p_i as "odds" w_i, then
# K(t) = N log ( sum_i w_i e^{t_i} ) - N log ( sum_j w_j )
#      = N log ( sum_i w_i (1 + e^{t_i} - 1) ) - N log ( sum_j w_j )
#      = N log ( sum_j w_j + sum_i w_i (e^{t_i} - 1) ) - N log ( sum_j w_j )
#      = N log (1 + sum_i w_i (e^{t_i} - 1) / sum_j w_j )
#      = N log1p ( sum_i w_i (e^{t_i} - 1) / sum_j w_j )
#      = N log1p ( sum_i w_i (expm1(t_i)) / sum_j w_j )

###      = N log1p ( sum_i w_i (exp(t_i) - 1) / sum_j w_j )
# NOTE: expm1(t_i) is numerically nicer, but we use exp(t_i)-1 for RTMB compatibility.


# NOTE: This class is single-block; replication across iid blocks is handled by iidReplicatesCGF()
# using block_size(parameter_vector) = length(parameter_vector) - 1.

MultinomialFamilyCGF <- R6::R6Class(
  classname = "MultinomialFamilyCGF",
  inherit   = CGF,

  private = list(

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    # K(t) via zm1 = exp(t) - 1 (RTMB-friendly)
    K_z1p = function(zm1, N_val, odds_val, odds_sum) {
      frac <- sum(odds_val * zm1) / odds_sum
      N_val * log1p(frac)
    },

    # v(t) = softmax in odds-parameterization:
    #   v_i = odds_i * exp(t_i) / sum_j odds_j * exp(t_j)
    v_from_t = function(tblock, odds_val) {
      numer <- odds_val * exp(tblock)
      numer / sum(numer)
    },

    .check_block = function(tvec, parameter_vector, where = "") {
      d <- length(parameter_vector) - 1
      if (d < 1L) {
        stop("MultinomialFamilyCGF: 'parameter_vector' must have length >= 2 (c(N, odds[1:d])).")
      }
      if (length(tvec) != d) {
        msg <- paste0(
          "MultinomialFamilyCGF", if (nzchar(where)) paste0("::", where) else "",
          ": length(tvec) must equal d = length(parameter_vector) - 1. ",
          "Got length(tvec) = ", length(tvec), ", d = ", d, "."
        )
        stop(msg)
      }
      invisible(d)
    },

    .check_Q = function(Q, d, where = "") {
      if (is.null(dim(Q)) || length(dim(Q)) != 2L) {
        stop("MultinomialFamilyCGF", if (nzchar(where)) paste0("::", where) else "",
             ": expected a matrix-like Q with dim().")
      }
      if (nrow(Q) != d || ncol(Q) != d) {
        stop("MultinomialFamilyCGF", if (nzchar(where)) paste0("::", where) else "",
             ": dimension mismatch: expected Q to be ", d, "x", d,
             ", got ", nrow(Q), "x", ncol(Q), ".")
      }
      invisible(TRUE)
    },

    # ------------------------------------------------------------------
    # Default methods (single block)
    # ------------------------------------------------------------------

    # parameter_vector = c(N, odds_1, ..., odds_d)
    # p_i = odds_i / sum(odds)
    # d = length(parameter_vector) - 1
    ###
    # K(t) = N * log( sum_i p_i exp(t_i) )
    # where p_i = odds_i / sum(odds)
    K_func_default = function(tvec, parameter_vector) {
      private$.check_block(tvec, parameter_vector, where = "K")
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      odds_sum <- sum(odds_val)
      zm1 <- exp(tvec) - 1
      private$K_z1p(zm1, N_val, odds_val, odds_sum)
    },

    # K'(t) = N * v(t)
    K1_func_default = function(tvec, parameter_vector) {
      private$.check_block(tvec, parameter_vector, where = "K1")
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)
      N_val * v
    },

    # K''(t) = N * (diag(v) - v v^T)
    K2_func_default = function(tvec, parameter_vector) {
      d <- private$.check_block(tvec, parameter_vector, where = "K2")
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      # Dense (RTMB-friendly). Note: this is singular for the full multinomial.
      N_val * (diag(v, nrow = d, ncol = d) - outer(v, v))
    },

    # 3rd-order cumulant tensor contraction:
    #   K3(u1,u2,u3) = sum_{i} v_i u1_i u2_i u3_i
    #                 - mu1 * sum_i v_i u2_i u3_i
    #                 - mu2 * sum_i v_i u1_i u3_i
    #                 - mu3 * sum_i v_i u1_i u2_i
    #                 + 2 mu1 mu2 mu3
    # then multiplied by N.
    K3operator_func_default = function(tvec, parameter_vector, u1, u2, u3) {
      d <- private$.check_block(tvec, parameter_vector, where = "K3operator")
      if (length(u1) != d || length(u2) != d || length(u3) != d) {
        stop("MultinomialFamilyCGF::K3operator: u1/u2/u3 must have length d.")
      }

      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      vu1 <- v * u1
      vu1s <- sum(vu1)
      vu2s <- sum(v * u2)
      vu3s <- sum(v * u3)

      u2u3 <- u2 * u3
      N_val * (sum(vu1 * u2u3) -
                 vu3s * sum(vu1 * u2) -
                 vu2s * sum(vu1 * u3) -
                 vu1s * sum(v * u2u3) +
                 2 * vu1s * vu2s * vu3s)
    },

    # 4th-order cumulant tensor contraction for u1,u2,u3,u4, multiplied by N.
    K4operator_func_default = function(tvec, parameter_vector, u1, u2, u3, u4) {
      d <- private$.check_block(tvec, parameter_vector, where = "K4operator")
      if (length(u1) != d || length(u2) != d || length(u3) != d || length(u4) != d) {
        stop("MultinomialFamilyCGF::K4operator: u1/u2/u3/u4 must have length d.")
      }

      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      vu1 <- v * u1
      vu2 <- v * u2
      vu3 <- v * u3
      vu4 <- v * u4

      vu1s <- sum(vu1)
      vu2s <- sum(vu2)
      vu3s <- sum(vu3)
      vu4s <- sum(vu4)

      u12 <- u1 * u2
      u34 <- u3 * u4

      vu12s <- sum(vu1 * u2) # = sum(v * u1 * u2)
      vu13s <- sum(vu1 * u3)
      vu14s <- sum(vu1 * u4)
      vu23s <- sum(vu2 * u3)
      vu24s <- sum(vu2 * u4)
      vu34s <- sum(vu3 * u4)

      vu123 <- u12 * vu3 # = v * u1 * u2 * u3

      N_val * (sum(vu123 * u4) -
                 vu4s * sum(vu123) -
                 vu3s * sum(u12 * vu4) -
                 vu2s * sum(u34 * vu1) -
                 vu1s * sum(u34 * vu2) -
                 vu12s * vu34s -
                 vu13s * vu24s -
                 vu14s * vu23s +
                 2 * (vu12s * vu3s * vu4s +
                        vu13s * vu2s * vu4s +
                        vu14s * vu2s * vu3s +
                        vu23s * vu1s * vu4s +
                        vu24s * vu1s * vu3s +
                        vu34s * vu1s * vu2s) -
                 6 * vu1s * vu2s * vu3s * vu4s)
    },

    # K4operatorAABB(t, Q1, Q2) for the multinomial:
    # This implementation assumes Q1 == Q2 (and ignores Q2).
    K4operatorAABB_func_default = function(tvec, parameter_vector, Q1, Q2) {
      d <- private$.check_block(tvec, parameter_vector, where = "K4operatorAABB")
      private$.check_Q(Q1, d, where = "K4operatorAABB")
      # Q2 is ignored by design (the calling code uses Q1 == Q2)
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      Qv <- Q1 %*% v
      vQv <- sum(v * Qv)

      res_double_indices <- sum(outer(v, v) * Q1 * Q1)

      diag_Q <- diag(Q1)
      tmp <- sum(v * diag_Q)

      N_val * (-2 * res_double_indices +
                 sum(v * diag_Q * (diag_Q - 2 * Qv - tmp + 2 * vQv)) -
                 sum(2 * v * Qv * diag_Q) +
                 sum(8 * v * Qv * Qv) +
                 2 * vQv * tmp -
                 6 * vQv * vQv)
    },

    # K3K3operatorAABBCC(t, Q1, Q2, Q3) for the multinomial.
    # This implementation assumes Q1 == Q2 == Q3 (and ignores Q2, Q3).
    K3K3operatorAABBCC_func_default = function(tvec, parameter_vector, Q1, Q2, Q3) {
      d <- private$.check_block(tvec, parameter_vector, where = "K3K3operatorAABBCC")
      private$.check_Q(Q1, d, where = "K3K3operatorAABBCC")
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      Qv  <- Q1 %*% v
      vQv <- sum(v * Qv)

      diag_Q <- diag(Q1)
      a <- v * diag_Q
      dvec <- v * Qv

      Q_a <- Q1 %*% a
      Q_d <- Q1 %*% dvec

      S1 <- sum(a    * Q_a)
      S2 <- 2 * sum(dvec * Q_a)
      S3 <- 2 * sum(a    * Q_d)
      S4 <- 4 * sum(dvec * Q_d)

      res_double_indices <- S1 - S2 - S3 + S4

      sum_v_diagQ_Qv <- sum(v * diag_Q * Qv)
      sum_v_diagQ    <- sum(v * diag_Q)
      sum_v_Qv_Qv    <- sum(v * Qv * Qv)

      N_val^2 * (res_double_indices +
                   sum_v_diagQ_Qv * (-sum_v_diagQ + 2 * vQv) +
                   sum_v_Qv_Qv    * ( 2 * sum_v_diagQ - 4 * vQv) +
                   sum_v_diagQ    * (-sum_v_diagQ_Qv +
                                       vQv * sum_v_diagQ +
                                       2 * sum_v_Qv_Qv -
                                       2 * vQv^2) +
                   2 * vQv        * (sum_v_diagQ_Qv -
                                       vQv * sum_v_diagQ -
                                       2 * sum_v_Qv_Qv +
                                       2 * vQv^2))
    },

    # K3K3operatorABCABC(t, Q1, Q2, Q3) for the multinomial.
    # This implementation assumes Q1 == Q2 == Q3 (and ignores Q2, Q3).
    K3K3operatorABCABC_func_default = function(tvec, parameter_vector, Q1, Q2, Q3) {
      d <- private$.check_block(tvec, parameter_vector, where = "K3K3operatorABCABC")
      private$.check_Q(Q1, d, where = "K3K3operatorABCABC")
      N_val    <- parameter_vector[1]
      odds_val <- parameter_vector[-1]
      v <- private$v_from_t(tvec, odds_val)

      Qv  <- Q1 %*% v
      vQv <- sum(v * Qv)

      len_v <- length(v)
      Qv_col <- matrix(Qv, nrow = len_v, ncol = len_v, byrow = FALSE)
      Qv_row <- matrix(Qv, nrow = len_v, ncol = len_v, byrow = TRUE)

      expression_matrix <- Q1^3 -
        3 * Q1^2 * Qv_col -
        3 * Q1^2 * Qv_row +
        3 * Q1^2 * vQv +
        6 * Q1 * Qv_col * Qv_row

      res_double_indices <- sum(outer(v, v) * expression_matrix)

      N_val^2 * (res_double_indices +
                   4  * sum(v * Qv^3) -
                   12 * vQv * sum(v * Qv^2) +
                   4  * vQv^3)
    },

    # Default correction term T(t) for SPA:
    func_Tfunc_default = function(tvec, parameter_vector) {
      Q <- solve(private$K2_func_default(tvec, parameter_vector))
      K3K3operatorABCABC_val <- private$K3K3operatorABCABC_func_default(tvec, parameter_vector, Q, Q, Q)
      K3K3operatorAABBCC_val <- private$K3K3operatorAABBCC_func_default(tvec, parameter_vector, Q, Q, Q)
      K4operatorAABB_val     <- private$K4operatorAABB_func_default(tvec, parameter_vector, Q, Q)
      K4operatorAABB_val / 8 - K3K3operatorAABBCC_val / 8 - K3K3operatorABCABC_val / 12
    },

    simulate_func_default = function(n, vector_length, parameter_vector, tvec = NULL, ...) {
      if (length(parameter_vector) < 2) {
        stop("MultinomialFamilyCGF$rsim: 'parameter_vector' must be c(N, odds[1:d]) with length >= 2.", call. = FALSE)
      }

      N_val <- parameter_vector[1]
      odds  <- parameter_vector[-1]
      d     <- length(odds)


      if (vector_length != d) {
        stop("MultinomialFamilyCGF$rsim: require vector_length == d = length(parameter_vector) - 1. ",
             "Got vector_length=", vector_length, ", d=", d, ".", call. = FALSE)
      }


      # validate N (multinomial size)
      if (!is.finite(N_val) || N_val < 0) stop("MultinomialFamilyCGF$rsim: N must be finite and >= 0.", call. = FALSE)

      tol <- 1e-8
      N_round <- round(N_val)
      if (abs(N_val - N_round) > tol) stop("MultinomialFamilyCGF$rsim: N must be (near) integer for simulation. Got N=", N_val, ".", call. = FALSE)

      N_int <- as.integer(N_round)


      if (any(!is.finite(odds)) || any(odds < 0)) stop("MultinomialFamilyCGF$rsim: odds/probabilities must be finite and >= 0.", call. = FALSE)

      odds_sum <- sum(odds)
      if (!is.finite(odds_sum) || odds_sum <= 0) {
        stop("MultinomialFamilyCGF$rsim: odds/probabilities must have a positive finite sum.", call. = FALSE)
      }

      # compute tilted probabilities if tvec is supplied:
      #   p_tilt ~ odds * exp(tvec)
      if (!is.null(tvec)) {
        if (length(tvec) != d) {
          stop("MultinomialFamilyCGF$rsim: if supplied, 'tvec' must have length d.", call. = FALSE)
        }
        if (any(!is.finite(tvec))) {
          stop("MultinomialFamilyCGF$rsim: 'tvec' must be finite.", call. = FALSE)
        }
        tmax <- max(tvec)
        w <- odds * exp(tvec - tmax)
      } else {
        w <- odds
      }

      wsum <- sum(w)

      p <- w / wsum

      # stats::rmultinom returns an integer matrix d x n
      stats::rmultinom(n = n, size = N_int, prob = p)
    }


  ),

  public = list(
    initialize = function(
      op_name = "MultinomialFamilyCGF",
      iidReps = "any",
      K = NULL,
      K1 = NULL,
      K2 = NULL,
      K3operator = NULL,
      K4operator = NULL,
      K2operator = NULL,
      K2operatorAK2AT = NULL,
      K4operatorAABB = NULL,
      K3K3operatorAABBCC = NULL,
      K3K3operatorABCABC = NULL,
      K4operatorAABB_factored = NULL,
      K3K3operatorAABBCC_factored = NULL,
      K3K3operatorABCABC_factored = NULL,
      ineq_constraint = NULL,
      analytic_tvec_hat = NULL,
      tilting_exponent = NULL,
      neg_ll = NULL,
      func_T = NULL,
      rsim = NULL,
      ...) {

      .check_iidReps(iidReps)
      if (!identical(iidReps, "any") && !(is.numeric(iidReps) && iidReps == 1L)) {
        stop("MultinomialFamilyCGF is single-block only. For replication use createMultinomialFamilyCGF(iidReps=...) or iidReplicatesCGF().")
      }

      # Defaults for the multinomial family
      final_K  <- if (is.null(K))  private$K_func_default  else K
      final_K1 <- if (is.null(K1)) private$K1_func_default else K1
      final_K2 <- if (is.null(K2)) private$K2_func_default else K2

      final_K3operator <- if (is.null(K3operator)) private$K3operator_func_default else K3operator
      final_K4operator <- if (is.null(K4operator)) private$K4operator_func_default else K4operator

      final_K4operatorAABB <- if (is.null(K4operatorAABB)) private$K4operatorAABB_func_default else K4operatorAABB
      final_K3K3operatorAABBCC <- if (is.null(K3K3operatorAABBCC)) private$K3K3operatorAABBCC_func_default else K3K3operatorAABBCC
      final_K3K3operatorABCABC <- if (is.null(K3K3operatorABCABC)) private$K3K3operatorABCABC_func_default else K3K3operatorABCABC

      # func_T_func <- if (is.null(func_Tfunc)) private$func_Tfunc_default else func_Tfunc
      final_func_T <- if(is.null(func_T)) private$func_Tfunc_default else func_T

      final_rsim <- if (is.null(rsim)) private$simulate_func_default else rsim




      super$initialize(
        K  = final_K,
        K1 = final_K1,
        K2 = final_K2,
        K3operator = final_K3operator,
        K4operator = final_K4operator,
        analytic_tvec_hat = analytic_tvec_hat,
        rsim = final_rsim,
        op_name = op_name,

        ineq_constraint = ineq_constraint,
        tilting_exponent = tilting_exponent,
        neg_ll  = neg_ll,
        func_T  = final_func_T,

        K2operator      = K2operator,
        K2operatorAK2AT = K2operatorAK2AT,

        K4operatorAABB      = final_K4operatorAABB,
        K3K3operatorAABBCC  = final_K3K3operatorAABBCC,
        K3K3operatorABCABC  = final_K3K3operatorABCABC,

        K4operatorAABB_factored     = K4operatorAABB_factored,
        K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
        K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,

        ...
      )


    }
  )
)


#' @keywords internal
createMultinomialFamilyCGF <- function(iidReps = "any",
                                       op_name = "MultinomialFamilyCGF",
                                       K = NULL,
                                       K1 = NULL,
                                       K2 = NULL,
                                       K3operator = NULL,
                                       K4operator = NULL,
                                       K2operator = NULL,
                                       K2operatorAK2AT = NULL,
                                       K4operatorAABB = NULL,
                                       K3K3operatorAABBCC = NULL,
                                       K3K3operatorABCABC = NULL,
                                       K4operatorAABB_factored = NULL,
                                       K3K3operatorAABBCC_factored = NULL,
                                       K3K3operatorABCABC_factored = NULL,
                                       ineq_constraint = NULL,
                                       analytic_tvec_hat = NULL,
                                       tilting_exponent = NULL,
                                       neg_ll = NULL,
                                       func_T = NULL,
                                       ...) {

  .check_iidReps(iidReps)

  base <- MultinomialFamilyCGF$new(
    op_name = op_name,
    K = K,
    K1 = K1,
    K2 = K2,
    K3operator = K3operator,
    K4operator = K4operator,
    K2operator = K2operator,
    K2operatorAK2AT = K2operatorAK2AT,
    K4operatorAABB = K4operatorAABB,
    K3K3operatorAABBCC = K3K3operatorAABBCC,
    K3K3operatorABCABC = K3K3operatorABCABC,
    K4operatorAABB_factored = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    ineq_constraint = ineq_constraint,
    analytic_tvec_hat = analytic_tvec_hat,
    tilting_exponent = tilting_exponent,
    neg_ll = neg_ll,
    func_T = func_T,
    ...
  )

  # Automatic block size:
  #   d(param) = length(param) - 1
  bs_fun <- function(parameter_vector) as.integer(length(parameter_vector) - 1L)
  attr(bs_fun, "label") <- "length(param_vector)-1"

  iidReplicatesCGF(cgf = base, iidReps = iidReps, block_size = bs_fun)
}
