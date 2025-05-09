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





#' MultinomialFamilyCGF Class
#'
#' @description
#' A subclass of `CGF` for the multinomial distribution. It allows the user to
#' optionally provide custom methods for `K, K1, K2, K3operator, K4operator` (the
#' five required methods), and also for any of the optional methods like
#' `tilting_exponent`, `neg_ll`, etc. If the user does not supply a function
#' (i.e., passes `NULL`), then a default multinomial-specific implementation is
#' used (for the five required methods). For any optional methods not supplied,
#' we defer to the base class defaults.
#'
#' @details
#' The parameter vector is assumed to be \eqn{(N, \omega_1, \ldots, \omega_d)},
#' where \eqn{N} is the total count, and \eqn{\omega_i} are "odds" that normalize
#' to probabilities \eqn{p_i = \omega_i / \sum_j \omega_j}.
#'
#' We implement:
#' \itemize{
#'   \item K(t) = N * log1p( sum(\omega_i [exp(t_i)-1]) / sum(\omega_i) ).
#'   \item K1(t), K2(t), K3operator(t), K4operator(t) 
#'         all derived from standard multinomial expansions.
#'   \item \code{func_T} (optional) for the first-order saddlepoint correction,
#'         referencing the built-in private factorized operators if you wish.
#' }
#'
#'
#' @noRd
MultinomialFamilyCGF <- R6::R6Class(
  classname = "MultinomialFamilyCGF",
  inherit   = CGF,   
  
  private = list(
    .iidReps = NULL,
    check_tvec_length = function(tvec, param) {
      if (!is.null(private$.iidReps)) {
        d <- length(param) - 1
        if (length(tvec) != d * private$.iidReps) {
          stop(sprintf(
            "MultinomialFamilyCGF: length(tvec)=%d != d(%d)*iidReps(%d).",
            length(tvec), d, private$.iidReps
          ))
        }
      }
    },
    
    # ------------------------------------------------------------------
    # Helper Functions
    # 
    # We'll store stable expansions for K, K1, K2, etc., 
    # possibly using expm1() and log1p() for numeric stability.
    # 
    # In the "param" vector, we will interpret param[1] = N,
    # and param[2..end] = odds vector (w_i).
    # ------------------------------------------------------------------
    
    
    #    We interpret 'odds_val' as w_i, 
    #    p_i = w_i / sum(w_j).
    #    Then sum_i [ w_i ( e^{t_i} - 1 ) ] / sum_j w_j 
    #    is the fraction inside log1p(...).
    K_z1p = function(zm1, N_val, odds_val, odds_sum) {
      # sum(odds_val * zm1) => sum_i [ w_i (e^{t_i}-1) ]
      # then divide by sum(w_i) => fraction => log1p
      frac <- sum(odds_val * zm1) / odds_sum
      N_val * log1p(frac)
    },
    
    # The tilted distribution v_i = p_i e^{t_i} / sum_j p_j e^{t_j}.
    # If p_i = w_i / sum(w_i), we do: 
    #    numerator = w_i * exp(t_i), 
    #    denom = sum_j [w_j * exp(t_j)].
    v_from_t = function(tvec, odds_val) {
      tmp <- odds_val * exp(tvec)
      tmp / sum(tmp)
    },
    
    
    #-------------------------------------------------------------------
    # Default versions of K, K1, K2, K3operator, K4operator
    # (the 5 required methods). Users can override by passing non-NULL.
    #-------------------------------------------------------------------
    K_func_default = function(tvec, param) {
      # private$check_tvec_length(tvec, param) 
      d <- length(param) - 1  # dimension for one block
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      odds_sum <- sum(odds_val)
      
      # Reshape tvec into a matrix with 'd' rows and 'nblocks' columns
      tmat <- matrix(tvec, nrow = d, ncol = nblocks)
      
      # Apply function to each column (block)
      block_values <- apply(tmat, 2, function(tblock) {
        zm1 <- exp(tblock) - 1
        private$K_z1p(zm1, N_val, odds_val, odds_sum)
      })
      
      sum(block_values)
    },
    

    # K1(tvec, param) = gradient wrt tvec
    # = N * v_from_t(tvec, p). 
    # v must be the normalised vector of probabilities for the tilted distribution
    # If p_i = w_i / sum(w_j), then v_i = w_i e^{t_i}/ sum_j [w_j e^{t_j}].
    K1_func_default = function(tvec, param) {
      # private$check_tvec_length(tvec, param) 
      d <- length(param) - 1  # dimension for a single block
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      tmat <- matrix(tvec, nrow = d, ncol = nblocks)
      
      block_results <- apply(tmat, 2, function(tblock) {
        v <- private$v_from_t(tblock, odds_val)
        N_val * v
      })
      
      as.vector(block_results)
    },
    
    # K2(tvec, param) = N [ diag(v) - v v^T ]
    # That is a d x d matrix.
    K2_func_default = function(tvec, param) {
      # private$check_tvec_length(tvec, param) 
      d <- length(param) - 1
      N_val    <- param[1]
      odds_val <- param[-1]
      
      n <- length(tvec)
      nblocks <- n / d
      
      big_mat <- Matrix::Matrix(0, nrow = n, ncol = n) * param[1]  ##### possibly sparse???
      
      for (i in seq_len(nblocks)) {
        idx_start <- (i - 1) * d + 1
        idx_end   <- i * d
        
        # Extract the portion of tvec belonging to block i
        tblock <- tvec[idx_start:idx_end]
        
        # Compute the block-level Hessian
        v    <- private$v_from_t(tblock, odds_val)
        vvT  <- outer(v, v)
        H_i  <- N_val * (diag(v, nrow = length(v)) - vvT)  # d x d
        
        big_mat[idx_start:idx_end, idx_start:idx_end] <- H_i
      }
      big_mat
    },
    
    # K3operator(tvec, u1, u2, u3, param).
    # The 3rd derivative "operator" at tvec, applied to vectors u1, u2, u3.
    K3operator_func_default = function(tvec, param, u1, u2, u3) {
      # private$check_tvec_length(tvec, param) 
      d <- length(param) - 1
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      
      # Reshape tvec and the vectors u1, u2, u3
      tmat  <- matrix(tvec, nrow = d, ncol = nblocks)
      u1mat <- matrix(u1,   nrow = d, ncol = nblocks)
      u2mat <- matrix(u2,   nrow = d, ncol = nblocks)
      u3mat <- matrix(u3,   nrow = d, ncol = nblocks)
      
      block_vals <- apply(seq_len(nblocks), 1, function(bi) {
        tblock <- tmat[, bi]
        u1b    <- u1mat[, bi]
        u2b    <- u2mat[, bi]
        u3b    <- u3mat[, bi]
        
        v <- private$v_from_t(tblock, odds_val)
        
        vu1  <- v * u1b
        vu1s <- sum(vu1)
        vu2s <- sum(v * u2b)
        vu3s <- sum(v * u3b)
        u2u3 <- u2b * u3b
        
        N_val * (
          sum(vu1 * u2u3) 
          - vu3s * sum(vu1 * u2b) 
          - vu2s * sum(vu1 * u3b) 
          - vu1s * sum(v * u2u3)
          + 2 * vu1s * vu2s * vu3s
        )
      })
      
      sum(block_vals)
    },
    
    K4operator_func_default = function(tvec, param, u1, u2, u3, u4) {
      # private$check_tvec_length(tvec, param) 
      d <- length(param) - 1
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      
      tmat  <- matrix(tvec, nrow = d, ncol = nblocks, byrow = FALSE)
      u1mat <- matrix(u1,   nrow = d, ncol = nblocks, byrow = FALSE)
      u2mat <- matrix(u2,   nrow = d, ncol = nblocks, byrow = FALSE)
      u3mat <- matrix(u3,   nrow = d, ncol = nblocks, byrow = FALSE)
      u4mat <- matrix(u4,   nrow = d, ncol = nblocks, byrow = FALSE)
      
      block_vals <- apply(seq_len(nblocks), 1, function(i) {
        tblock <- tmat[, i]
        u1b    <- u1mat[, i]
        u2b    <- u2mat[, i]
        u3b    <- u3mat[, i]
        u4b    <- u4mat[, i]
        
        v <- private$v_from_t(tblock, odds_val)
        
        vu1  <- v * u1b
        vu2  <- v * u2b
        vu3  <- v * u3b
        vu4  <- v * u4b
        vu1s <- sum(vu1)
        vu2s <- sum(vu2)
        vu3s <- sum(vu3)
        vu4s <- sum(vu4)
        u12  <- u1b * u2b
        u34  <- u3b * u4b
        
        vu12s <- sum(vu1 * u2b)
        vu13s <- sum(vu1 * u3b)
        vu14s <- sum(vu1 * u4b)
        vu23s <- sum(vu2 * u3b)
        vu24s <- sum(vu2 * u4b)
        vu34s <- sum(vu3 * u4b)
        vu123 <- u12 * vu3
        
        N_val*(
          sum(vu123*u4b) - vu4s*sum(vu123) - vu3s*sum(u12 * vu4) - vu2s*sum(u34 * vu1)
          - vu1s*sum(u34 * vu2) - vu12s*vu34s - vu13s*vu24s - vu14s*vu23s
          + 2*(vu12s*vu3s*vu4s + vu13s*vu2s*vu4s + vu14s*vu2s*vu3s
               + vu23s*vu1s*vu4s + vu24s*vu1s*vu3s + vu34s*vu1s*vu2s)
          - 6*vu1s*vu2s*vu3s*vu4s
        )
        
      })
      
      sum(block_vals)
      },
    
    K4operatorAABB_func_default = function(tvec, param, Q1, Q2) {
      # private$check_tvec_length(tvec, param) 
      # # v <- private$v_from_t(tvec, odds_val)
      # 
      # # Q1v  <- Q1 %*% v
      # # Q2v  <- Q2 %*% v
      # # vQ1v <- sum(v * Q1v)
      # # vQ2v <- sum(v * Q2v)
      # # res_double_indices <- sum(outer(v, v) * Q1 * Q2)
      # # tmp <- sum(v * diag(Q2))
      # # N_val * (
      # #   -2*res_double_indices +
      # #     sum(v * diag(Q1) * (diag(Q2) - 2*Q2v - tmp + 2*vQ2v)) -
      # #     sum(2*v*Q1v*diag(Q2)) +
      # #     sum(8*v*Q1v*Q2v) +
      # #     2*vQ1v*tmp -
      # #     6*vQ1v*vQ2v
      # # )
      
      d <- length(param) - 1
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      
      process_block <- function(i) {
        idx_start <- (i - 1) * d + 1
        idx_end   <- i * d
        
        tblock <- tvec[idx_start:idx_end]
        v <- private$v_from_t(tblock, odds_val)
        Q1block <- Q1[idx_start:idx_end, idx_start:idx_end]
        
        # Q1 == Q2
        diag_Q1 <- diag(Q1block)
        Q1v  <- Q1block %*% v  
        vQ1v <- sum(v * Q1v)
        
        res_double_indices <- sum(outer(v, v) * Q1block * Q1block)
        tmp <- sum(v * diag_Q1)
        
        N_val * (
          -2*res_double_indices +
            sum(v*diag_Q1 *(diag_Q1 - 2*Q1v - tmp + 2*vQ1v)) -
            sum(2*v*Q1v*diag_Q1) +
            sum(8*v*Q1v*Q1v) +
            2*vQ1v*tmp -
            6*vQ1v*vQ1v
        )

      }
      operator_values <- sapply(seq_len(nblocks), process_block)
      sum(operator_values)
      
    },
    
    K3K3operatorAABBCC_func_default = function(tvec, param, Q1, Q2, Q3) {
      # private$check_tvec_length(tvec, param) 
      
      # Q1v <- Q1 %*% v
      # Q2v <- Q2 %*% v
      # Q3v <- Q3 %*% v
      # vQ1v <- sum(v * Q1v)
      # vQ2v <- sum(v * Q2v)
      # vQ3v <- sum(v * Q3v)
      # diag_Q1 <- diag(Q1)
      # diag_Q3 <- diag(Q3)
      # a <- v * diag_Q1
      # b <- v * diag_Q3
      # d <- v * Q1v
      # c <- v * Q3v
      # Q2_b <- Q2 %*% b          # Q2 %*% (v * diag_Q3)
      # Q2_c <- Q2 %*% c          # Q2 %*% (v * Q3v)
      # S1 <- sum(a * Q2_b)
      # S2 <- 2 * sum(d * Q2_b)
      # S3 <- 2 * sum(a * Q2_c)
      # S4 <- 4 * sum(d * Q2_c)
      # res_double_indices <- S1 - S2 - S3 + S4
      # sum_v_diagQ3_Q2v   <- sum(v * diag_Q3 * Q2v)
      # sum_v_diagQ1       <- sum(v * diag_Q1)
      # sum_v_Q2v_Q3v      <- sum(v * Q2v * Q3v)
      # sum_v_diagQ3       <- sum(v * diag_Q3)
      # sum_v_diagQ1_Q2v   <- sum(v * diag_Q1 * Q2v)
      # sum_v_Q1v_Q2v      <- sum(v * Q1v * Q2v)
      # 
      # N_val^2 * (
      #   res_double_indices +
      #     sum_v_diagQ3_Q2v * (-sum_v_diagQ1 + 2 * vQ1v) +
      #     sum_v_Q2v_Q3v * (2 * sum_v_diagQ1 - 4 * vQ1v) +
      #     sum_v_diagQ3 * (
      #       -sum_v_diagQ1_Q2v +
      #         vQ2v * sum_v_diagQ1 +
      #         2 * sum_v_Q1v_Q2v -
      #         2 * vQ1v * vQ2v
      #     ) +
      #     2 * vQ3v * (
      #       sum_v_diagQ1_Q2v -
      #         vQ2v * sum_v_diagQ1 -
      #         2 * sum_v_Q1v_Q2v +
      #         2 * vQ1v * vQ2v
      #     )
      # )
      
      
      
      # Q1 == Q2 == Q3
      Q <- Q1
      d <- length(param) - 1
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      
      process_block <- function(i) {
        idx_start <- (i - 1) * d + 1
        idx_end   <- i * d
        
        tblock <- tvec[idx_start:idx_end]
        v <- private$v_from_t(tblock, odds_val)
        Qblock <- Q[idx_start:idx_end, idx_start:idx_end]
        
        Qv <- Qblock %*% v
        vQv <- sum(v * Qv)
        diag_Q <- diag(Qblock)
        a <- v * diag_Q
        d <- v * Qv
        Q_a <- Qblock %*% a          # Equivalent to Qblock %*% (v * diag(Qblock))
        Q_d <- Qblock %*% d          # Equivalent to Qblock %*% (v * Qv)
        S1 <- sum(a * Q_a)
        S2 <- 2 * sum(d * Q_a)
        S3 <- 2 * sum(a * Q_d)
        S4 <- 4 * sum(d * Q_d)
        res_double_indices <- S1 - S2 - S3 + S4
        sum_v_diagQ_Qv  <- sum(v * diag_Q * Qv)    
        sum_v_diagQ     <- sum(v * diag_Q)        
        sum_v_Qv_Qv     <- sum(v * Qv * Qv)        
        
        N_val^2*(
          res_double_indices +
            sum_v_diagQ_Qv * (-sum_v_diagQ + 2 * vQv) +
            sum_v_Qv_Qv * (2 * sum_v_diagQ - 4 * vQv) +
            sum_v_diagQ * (
              -sum_v_diagQ_Qv +
                vQv * sum_v_diagQ +
                2 * sum_v_Qv_Qv -
                2 * vQv^2
            ) +
            2 * vQv * (
              sum_v_diagQ_Qv -
                vQv * sum_v_diagQ -
                2 * sum_v_Qv_Qv +
                2 * vQv^2
            )
        )

      }
      operator_values <- sapply(seq_len(nblocks), process_block)
      sum(operator_values)
    },
    
    K3K3operatorABCABC_func_default = function(tvec, param, Q1, Q2, Q3) {
      # private$check_tvec_length(tvec, param) 
      
      # Q1v <- as.vector(Q1 %*% v)
      # Q2v <- as.vector(Q2 %*% v)
      # Q3v <- as.vector(Q3 %*% v)
      # 
      # vQ1v <- sum(v * Q1v)
      # vQ2v <- sum(v * Q2v)
      # vQ3v <- sum(v * Q3v)
      # 
      # # res_double_indices = 0
      # # for(j1 in 1:length(v)){
      # #   for(j2 in 1:length(v)){
      # #     res_double_indices = res_double_indices + (v[j1]*v[j2]) * (Q1[j1,j2]*Q2[j1,j2] * (Q3[j1,j2] - Q3v[j1] - Q3v[j2] + vQ3v) +
      # #                                                                  Q2[j1,j2]*Q3[j1,j2] * (-Q1v[j1] - Q1v[j2] + vQ1v) +
      # #                                                                  Q2[j1,j2]*(Q1v[j1]*Q3v[j2] + Q1v[j2]*Q3v[j1] ) +
      # #                                                                  Q1[j1,j2]*Q3[j1,j2]*(-Q2v[j1] - Q2v[j2] ) +
      # #                                                                  Q1[j1,j2]*(Q2v[j1]*Q3v[j2] + Q2v[j2]*Q3v[j1] ) +
      # #                                                                  Q3[j1,j2]*(Q1v[j1]*Q2v[j2] + Q1v[j2]*Q2v[j1] + Q1[j1,j2]*vQ2v))
      # #   }
      # # }
      # 
      # len_v <- length(v)
      # Q1v_col <- Matrix(as.vector(Q1v), nrow = len_v, ncol = len_v, byrow = FALSE)
      # Q1v_row <- Matrix(Q1v, nrow = len_v, ncol = len_v, byrow = TRUE)
      # Q2v_col <- Matrix(Q2v, nrow = len_v, ncol = len_v, byrow = FALSE)
      # Q2v_row <- Matrix(Q2v, nrow = len_v, ncol = len_v, byrow = TRUE)
      # Q3v_col <- Matrix(Q3v, nrow = len_v, ncol = len_v, byrow = FALSE)
      # Q3v_row <- Matrix(Q3v, nrow = len_v, ncol = len_v, byrow = TRUE)
      # 
      # 
      # term1 <- Q1 * Q2 * (Q3 - Q3v_col - Q3v_row + vQ3v)
      # term2 <- Q2 * Q3 * (-Q1v_col - Q1v_row + vQ1v)
      # term3 <- Q2 * (Q1v_col * Q3v_row + Q1v_row * Q3v_col)
      # term4 <- Q1 * Q3 * (-Q2v_col - Q2v_row)
      # term5 <- Q1 * (Q2v_col * Q3v_row + Q2v_row * Q3v_col)
      # term6 <- Q3 * (Q1v_col * Q2v_row + Q1v_row * Q2v_col + Q1 * vQ2v)
      # 
      # 
      # expression_matrix <- term1 + term2 + term3 + term4 + term5 + term6
      # outer_vv <- outer(v, v)
      # res_double_indices <- sum(outer_vv * expression_matrix)
      # 
      # 
      # N_val^2 * (
      #   res_double_indices
      #   + 4*sum(v*Q1v*Q2v*Q3v)
      #   - 4*vQ3v*sum(v*Q1v*Q2v)
      #   - 4*vQ2v*sum(v*Q1v*Q3v)
      #   - 4*vQ1v*sum(v*Q2v*Q3v)
      #   + 4*vQ1v*vQ2v*vQ3v
      # )
      

      
      # Q1 == Q2 == Q3
      Q <- Q1
      d <- length(param) - 1
      nblocks <- length(tvec) / d
      
      N_val    <- param[1]
      odds_val <- param[-1]
      
      process_block <- function(i) {
        idx_start <- (i - 1) * d + 1
        idx_end   <- i * d
        
        tblock <- tvec[idx_start:idx_end]
        v <- private$v_from_t(tblock, odds_val)
        Qblock <- Q[idx_start:idx_end, idx_start:idx_end]
        
        
        Qv <- as.vector(Qblock %*% v) 
        vQv <- sum(v * Qv)               
        
        len_v <- length(v)
        Qv_col <- matrix(Qv, nrow = len_v, ncol = len_v, byrow = FALSE)
        Qv_row <- matrix(Qv, nrow = len_v, ncol = len_v, byrow = TRUE)
        
        expression_matrix <- Qblock^3 - 
                3 * Qblock^2 * Qv_col - 
                3 * Qblock^2 * Qv_row + 
                3 * Qblock^2 * vQv + 
                6 * Qblock * Qv_col * Qv_row
        
        outer_vv <- outer(v, v)
        res_double_indices <- sum(outer_vv * expression_matrix)
        
        N_val^2 * (
          res_double_indices +
            4 * sum(v * Qv^3) -
            12 * vQv * sum(v * Qv^2) +
            4 * (vQv)^3
        )
        
      }
      operator_values <- sapply(seq_len(nblocks), process_block)
      sum(operator_values)
    },
    
    #### We avoid the factored forms for now 
    func_Tfunc_default = function(tvec, param) {
      Q <- solve(private$K2_func_default(tvec, param))
      K3K3operatorABCABC_val <- private$K3K3operatorABCABC_func_default(tvec, param, Q, Q, Q)
      K3K3operatorAABBCC_val <- private$K3K3operatorAABBCC_func_default(tvec, param, Q, Q, Q)
      K4operatorAABB_val <- private$K4operatorAABB_func_default(tvec, param, Q, Q)
      K4operatorAABB_val/8 - K3K3operatorAABBCC_val/8 - K3K3operatorABCABC_val/12
    }
  ),
    
    
  
  public = list(
    
    # @description
    # Constructor. The user may supply any or all of the 5 required functions
    # (K, K1, K2, K3operator, K4operator). If any is `NULL`, the default
    # multinomial version is used. For the optional methods, passing them as
    # arguments will override the base defaults. If they're `NULL`, we keep
    # the base class's default.
    # 
    # @param K_func,K1_func,K2_func,K3operator_func,K4operator_func (Required methods)
    #  If `NULL`, use the multinomial default. If a function, override the default.
    # @param ineq_constraint_func Optional.
    # @param analytic_tvec_hat_func Optional.
    # @param tilting_exponent_func Optional.
    # @param neg_ll_func Optional.
    # @param func_T_func Optional.
    # @param K4operatorAABB_func,K3K3operatorAABBCC_func,K3K3operatorABCABC_func Optional.
    # @param K4operatorAABB_factored_func,K3K3operatorAABBCC_factored_func,K3K3operatorABCABC_factored_func Optional.
    # @param K2operator_func,K2operatorAK2AT_func Optional.
    # @param op_name Character string label.
    # @param ... Additional arguments for future use or to pass to the CGF's \code{initialize}.
    initialize = function(
      K_func = NULL,
      K1_func = NULL,
      K2_func = NULL,
      K3operator_func = NULL,
      K4operator_func = NULL,
      
      ineq_constraint_func = NULL,
      analytic_tvec_hat_func = NULL,
      tilting_exponent_func = NULL,
      neg_ll_func = NULL,
      func_T_func = NULL,
      
      K4operatorAABB_func = NULL,
      K3K3operatorAABBCC_func = NULL,
      K3K3operatorABCABC_func = NULL,
      
      K4operatorAABB_factored_func = NULL,
      K3K3operatorAABBCC_factored_func = NULL,
      K3K3operatorABCABC_factored_func = NULL,
      
      K2operator_func = NULL,
      K2operatorAK2AT_func = NULL,
      
      op_name = "UnnamedOperation",
      iidReps = "any",
      ...
    ) {
      
      if (is.character(iidReps) && length(iidReps)==1 && tolower(iidReps)=="any") {
        private$.iidReps <- NULL
      }else {
        if (!is.numeric(iidReps) || length(iidReps)!=1 || iidReps<1 || iidReps!=as.integer(iidReps) || is.infinite(iidReps)) {
          stop("'iidReps' must be 'any' or a positive integer.")
        }
        private$.iidReps <- as.integer(iidReps)
      }
      
      # If the user-supplied function is NULL, use the child's default. 
      # Otherwise use the user's function.
      
      final_K_func  <- if(is.null(K_func))  private$K_func_default  else K_func
      final_K1_func <- if(is.null(K1_func)) private$K1_func_default else K1_func
      final_K2_func <- if(is.null(K2_func)) private$K2_func_default else K2_func
      final_K3op_func <- if(is.null(K3operator_func)) private$K3operator_func_default else K3operator_func
      final_K4op_func <- if(is.null(K4operator_func)) private$K4operator_func_default else K4operator_func
      final_K4opAABB_func <- if(is.null(K4operatorAABB_func)) private$K4operatorAABB_func_default else K4operatorAABB_func
      final_K3K3opAABBCC_func <- if(is.null(K3K3operatorAABBCC_func)) private$K3K3operatorAABBCC_func_default else K3K3operatorAABBCC_func
      final_K3K3opABCABC_func <- if(is.null(K3K3operatorABCABC_func)) private$K3K3operatorABCABC_func_default else K3K3operatorABCABC_func
      final_func_T_func <- if(is.null(func_T_func)) private$func_Tfunc_default else func_T_func
      
      
      super$initialize(
        # The 5 required:
        K_func         = final_K_func,
        K1_func        = final_K1_func,
        K2_func        = final_K2_func,
        K3operator_func= final_K3op_func,
        K4operator_func= final_K4op_func,
        
        # Optional methods:
        ineq_constraint_func = ineq_constraint_func,
        analytic_tvec_hat_func = analytic_tvec_hat_func,
        tilting_exponent_func  = tilting_exponent_func,
        neg_ll_func            = neg_ll_func,
        func_T_func            = final_func_T_func,
        
        K4operatorAABB_func         = final_K4opAABB_func,
        K3K3operatorAABBCC_func     = final_K3K3opAABBCC_func,
        K3K3operatorABCABC_func     = final_K3K3opABCABC_func,
        
        K4operatorAABB_factored_func    = K4operatorAABB_factored_func,
        K3K3operatorAABBCC_factored_func= K3K3operatorAABBCC_factored_func,
        K3K3operatorABCABC_factored_func= K3K3operatorABCABC_factored_func,
        
        K2operator_func       = K2operator_func,
        K2operatorAK2AT_func  = K2operatorAK2AT_func,
        
        op_name = op_name,
        ...
      )
      
      
    }
  )
)






#' Create a MultinomialFamilyCGF object
#'
#' @description
#' This factory function instantiates a CGF object related to the multinomial distribution.
#' It accepts arguments for all required and optional CGF methods. 
#' If a required method (K_func, K1_func, etc.) is `NULL`, the class's built-in 
#' multinomial default is used. For optional methods that are `NULL`, we rely on 
#' the base CGF defaults.
#'
#' @param K_func,K1_func,K2_func,K3operator_func,K4operator_func The five required methods; if `NULL`, 
#'   use the multinomial default in the child class.
#' @param ineq_constraint_func Optional function (default `NULL` => no constraints).
#' @param analytic_tvec_hat_func Optional function.
#' @param tilting_exponent_func,neg_ll_func,func_T_func Optional overrides.
#' @param K4operatorAABB_func,K3K3operatorAABBCC_func,K3K3operatorABCABC_func Optional overrides.
#' @param K4operatorAABB_factored_func,K3K3operatorAABBCC_factored_func,K3K3operatorABCABC_factored_func Optional overrides.
#' @param K2operator_func,K2operatorAK2AT_func Optional overrides.
#' @param op_name Character string operation label. Default "UnnamedOperation".
#' @param iidReps Either \code{"any"} or a positive integer specifying how many i.i.d. blocks are expected. Defaults to \code{"any"}.
#' @param ... Additional arguments for future use or to pass to the base CGF's \code{initialize}.
#'
#' @return A CGF object.
#' @noRd
createMultinomialFamilyCGF <- function(
    K_func = NULL,
    K1_func = NULL,
    K2_func = NULL,
    K3operator_func = NULL,
    K4operator_func = NULL,
    
    ineq_constraint_func = NULL,
    analytic_tvec_hat_func = NULL,
    tilting_exponent_func = NULL,
    neg_ll_func = NULL,
    func_T_func = NULL,
    
    K4operatorAABB_func = NULL,
    K3K3operatorAABBCC_func = NULL,
    K3K3operatorABCABC_func = NULL,
    
    K4operatorAABB_factored_func = NULL,
    K3K3operatorAABBCC_factored_func = NULL,
    K3K3operatorABCABC_factored_func = NULL,
    
    K2operator_func = NULL,
    K2operatorAK2AT_func = NULL,
    iidReps = "any",
    
    op_name = "UnnamedOperation",
    ...
) {
  MultinomialFamilyCGF$new(
    K_func = K_func,
    K1_func = K1_func,
    K2_func = K2_func,
    K3operator_func = K3operator_func,
    K4operator_func = K4operator_func,
    
    ineq_constraint_func = ineq_constraint_func,
    analytic_tvec_hat_func = analytic_tvec_hat_func,
    tilting_exponent_func  = tilting_exponent_func,
    neg_ll_func            = neg_ll_func,
    func_T_func            = func_T_func,
    
    K4operatorAABB_func         = K4operatorAABB_func,
    K3K3operatorAABBCC_func     = K3K3operatorAABBCC_func,
    K3K3operatorABCABC_func     = K3K3operatorABCABC_func,
    
    K4operatorAABB_factored_func     = K4operatorAABB_factored_func,
    K3K3operatorAABBCC_factored_func = K3K3operatorAABBCC_factored_func,
    K3K3operatorABCABC_factored_func = K3K3operatorABCABC_factored_func,
    
    K2operator_func      = K2operator_func,
    K2operatorAK2AT_func = K2operatorAK2AT_func,
    
    op_name = op_name,
    iidReps = iidReps,
    ...
  )
}



























































