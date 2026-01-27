# --------------------------------------------------------------------
# File: A-CGF_base.R
#
# PURPOSE & INSTRUCTIONS:
#  1) This file defines a base CGF (cumulant generating function) class using R6.
#  2) There are five compulsory methods any CGF must provide:
#       - K(tvec, parameter_vector)
#       - K1(tvec, parameter_vector)
#       - K2(tvec, parameter_vector)
#       - K3operator(tvec, parameter_vector, v1, v2, v3)
#       - K4operator(tvec, parameter_vector, v1, v2, v3, v4)
#  3) The class also supports optional methods (tilting_exponent, neg_ll, func_T, etc.). Users
#     can supply them or rely on defaults.
#  4) Some methods are private (e.g., neg_ll, func_T, and operator "factored" forms). They remain
#     in the private environment, and are exposed via the public `.private_api` handle to avoid
#     cluttering the public method list.
#  6) The `createCGF()` factory function at the end makes it easy to instantiate a CGF without subclassing.
#
#  NOTE ON ENVIRONMENTS:
#  - All names you intend to overwrite in the constructor must be declared
#    in 'public = list(...)' or 'private = list(...)' to avoid "cannot add bindings
#    to a locked environment" errors.
#  - By default, this class is not aware of the "param_adaptor". If you need parameter adaptation, see adaptCGF().
#
#  EXTENSIBILITY:
#  - The 'initialize' method accepts '...' for additional named methods or overrides.
#  - End-users typically won't instantiate this class directly if they're using a specialized CGF,
#    but they can use createCGF(...) as a quick factory approach.
# --------------------------------------------------------------------


# ------------------------------------------------------------------------
# Helper function: check_fun_sig
# This check that the user-supplied function has the right signature
# Currently not used.
# ------------------------------------------------------------------------
check_fun_sig <- function(fn, expected_args) {
  if (!is.function(fn)) stop("Supplied object is not a function.")
  formals_list <- formals(fn)
  actual_args <- names(formals_list)
  if (!identical(actual_args, expected_args)) {
    stop(
      "Invalid function signature. Expected arguments: ",
      paste(expected_args, collapse = ", "),
      "; got: ",
      paste(actual_args, collapse = ", ")
    )
  }
}


#------------------------------------------------------------------------
# NOTE (log-determinants under RTMB):
# Earlier versions included an ADjoint-based helper to customize reverse-mode
# behaviour for log|det(K2)|. We rely on determinant() instead.
#
# Recommended pattern:
#   * use cgf$logdetK2(t, theta) whenever you need log|det(K2(t,theta))|
#   * use cgf$K2_solve(t, theta, rhs) whenever you need (K2)^{-1} rhs
#
# These have defaults based on determinant(cgf$K2(...), logarithm=TRUE) and
# solve(cgf$K2(...), rhs), but wrapper CGFs (iidReplicates / concatenation /
# multinomial, etc.) can override them to avoid ever materialising a huge global K2.
#------------------------------------------------------------------------


CGF_public_defaults <- list(
  K2operator = function(tvec, parameter_vector, x, y) {
    K2_val <- self$K2(tvec, parameter_vector)
    as.vector(t(x) %*% (K2_val %*% y))
  }
  ,
  K2operatorAK2AT = function(tvec, parameter_vector, A) {
    K2_val <- self$K2(tvec, parameter_vector)
    A %*% K2_val %*% t(A)
  }
  ,
  # Default behaviour uses solve(K2, rhs), but can be overridden for structure/speed.
  K2_solve = function(tvec, parameter_vector, rhs) {
    K2_val <- self$K2(tvec, parameter_vector)
    solve(K2_val, rhs)
  }
  ,
  # Default uses determinant() for RTMB compatibility.
  logdetK2 = function(tvec, parameter_vector) {
    K2_val <- self$K2(tvec, parameter_vector)
    determinant(K2_val, logarithm = TRUE)$modulus
  }
  ,
  # Simulation wrapper. Always present, but errors unless has_rsim == TRUE.
  rsim = function(n, vector_length, parameter_vector, tvec = NULL, flatten = FALSE, ...) {
    if (!isTRUE(self$has_rsim) || is.null(private$rsim_func)) {
      stop("This CGF does not implement simulation (no 'rsim' supplied).", call. = FALSE)
    }
    if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
      stop("'n' must be a positive integer.", call. = FALSE)
    }
    n <- as.integer(n)

    if (length(vector_length) != 1L || !is.finite(vector_length) || vector_length < 1L ||
        vector_length != as.integer(vector_length)) {
      stop("'vector_length' must be a positive integer.", call. = FALSE)
    }
    vector_length <- as.integer(vector_length)

    if (!is.logical(flatten) || length(flatten) != 1L || is.na(flatten)) {
      stop("'flatten' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(tvec)) {
      if (!is.numeric(tvec)) stop("'tvec' must be NULL or numeric.", call. = FALSE)
      if (length(tvec) != vector_length) {
        stop("'tvec' must have length == vector_length.", call. = FALSE)
      }
      if (any(!is.finite(tvec))) stop("'tvec' must be finite.", call. = FALSE)
    }

    out <- private$rsim_func(n, vector_length, parameter_vector, tvec, ...)

    # Enforce: either a numeric vector of length n*vector_length, or a (vector_length x n) numeric matrix
    if (is.null(dim(out))) {
      if (!is.numeric(out)) stop("rsim must return a numeric vector or matrix.", call. = FALSE)
      if (length(out) != n * vector_length) {
        stop(
          "rsim returned a vector of length ", length(out),
          ", expected ", n * vector_length, " (= n * vector_length).",
          call. = FALSE
        )
      }
      out <- matrix(out, nrow = vector_length, ncol = n)
    } else {
      out <- as.matrix(out)
      if (!is.numeric(out)) stop("rsim must return a numeric vector or matrix.", call. = FALSE)
      if (nrow(out) != vector_length || ncol(out) != n) {
        stop(
          "rsim returned a matrix with dim=", paste(dim(out), collapse = "x"),
          ", expected ", vector_length, "x", n, " (= vector_length x n).",
          call. = FALSE
        )
      }
    }

    if (flatten) return(as.numeric(out))
    out
  }
  ,
  K4operatorAABB = function(tvec, parameter_vector, Q1, Q2) {
    chol_Q1 <- chol(Q1)
    diag_Q1 <- diag(chol_Q1)
    d1 <- diag_Q1 * diag_Q1
    A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
    private$K4operatorAABB_factored(tvec, parameter_vector, A1, d1, A1, d1)
  }
  ,
  K3K3operatorAABBCC = function(tvec, parameter_vector, Q1, Q2, Q3) {
    chol_Q1 <- chol(Q1)
    diag_Q1 <- diag(chol_Q1)
    d1 <- diag_Q1 * diag_Q1
    A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
    private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1)
  }
  ,
  K3K3operatorABCABC = function(tvec, parameter_vector, Q1, Q2, Q3) {
    chol_Q1 <- chol(Q1)
    diag_Q1 <- diag(chol_Q1)
    d1 <- diag_Q1 * diag_Q1
    A1 <- t(chol_Q1) %*% diag(1/diag_Q1)
    private$K3K3operatorABCABC_factored(tvec, parameter_vector, A1, d1, A1, d1, A1, d1)
  }
  ,
  ineq_constraint = function(tvec, parameter_vector) {
    numeric(0)
  }
)


















CGF_private_defaults <- list(
  tilting_exponent = function(tvec, parameter_vector) {
    self$K(tvec, parameter_vector) - sum(tvec * self$K1(tvec, parameter_vector))
  }
  ,
  neg_ll = function(tvec, parameter_vector) {
    te <- private$tilting_exponent(tvec, parameter_vector)
    val_logdet <- self$logdetK2(tvec, parameter_vector)
    0.5 * val_logdet + 0.5 * length(tvec) * log(2*pi) - te
  }
  ,
  func_T = function(tvec, parameter_vector) {
    K2_val <- self$K2(tvec, parameter_vector)
    K2_inv <- solve(K2_val)
    chol_K2_inv <- chol(K2_inv)
    diag_K2_inv <- diag(chol_K2_inv)
    d <- diag_K2_inv * diag_K2_inv
    A <- t(chol_K2_inv) %*% diag(1/diag_K2_inv)

    K4_AABB   <- private$K4operatorAABB_factored(tvec, parameter_vector, A, d, A, d)
    K3K3_ABBC <- private$K3K3operatorAABBCC_factored(tvec, parameter_vector, A, d, A, d, A, d)
    K3K3_ABC  <- private$K3K3operatorABCABC_factored(tvec, parameter_vector, A, d, A, d, A, d)
    K4_AABB/8 - K3K3_ABBC/8 - K3K3_ABC/12
  }
  ,
  K4operatorAABB_factored = function(tvec, parameter_vector, A1, d1, A2, d2) {
    r1 <- length(d1)
    r2 <- length(d2)
    res <- 0
    for (m1 in seq_len(r1)) {
      for (m2 in seq_len(r2)) {
        res <- res + d1[m1]*d2[m2]*self$K4operator(
          tvec, parameter_vector, A1[,m1], A1[,m1], A2[,m2], A2[,m2]
        )
      }
    }
    res
  }
  ,
  K3K3operatorAABBCC_factored = function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
    r1 <- length(d1)
    r2 <- length(d2)
    r3 <- length(d3)
    res <- 0
    for (m2 in seq_len(r2)) {
      factor1 <- 0
      for (m1 in seq_len(r1)) {
        factor1 <- factor1 + d1[m1]*self$K3operator(tvec, parameter_vector, A1[,m1], A1[,m1], A2[,m2])
      }
      factor2 <- 0
      for (m3 in seq_len(r3)) {
        factor2 <- factor2 + d3[m3]*self$K3operator(tvec, parameter_vector, A2[,m2], A3[,m3], A3[,m3])
      }
      res <- res + d2[m2]*factor1*factor2
    }
    res
  }
  ,
  K3K3operatorABCABC_factored = function(tvec, parameter_vector, A1, d1, A2, d2, A3, d3) {
    r1 <- length(d1)
    r2 <- length(d2)
    r3 <- length(d3)
    message("The discrepancy option/compute.funcT has initiated a computation that may take a few moments...")
    res <- 0
    for (m1 in seq_len(r1)) {
      for (m2 in seq_len(r2)) {
        for (m3 in seq_len(r3)) {
          val <- self$K3operator(tvec, parameter_vector, A1[,m1], A2[,m2], A3[,m3])
          res <- res + d1[m1]*d2[m2]*d3[m3]*(val*val)
        }
      }
    }
    res
  }
)
















#' @noRd
CGF <- R6::R6Class(
  classname = "CGF",

  private = c(CGF_private_defaults, list(
    analytic_tvec_hat_func = NULL,
    rsim_func = NULL
  )),

  active = list(
    .private_api = function(value) {
      if (!missing(value)) stop("'.private_api' is read-only.", call. = FALSE)
      private
    }
  ),

  public = c(CGF_public_defaults, list(
    call_history = NULL,

    has_analytic_tvec_hat = FALSE,
    analytic_tvec_hat = NULL,

    has_rsim = FALSE,

    additional_methods = list(),

    # Required methods (set in initialize)
    K = function(tvec, parameter_vector) {
      stop("CGF$K is not initialized.", call. = FALSE)
    },
    K1 = function(tvec, parameter_vector) {
      stop("CGF$K1 is not initialized.", call. = FALSE)
    },
    K2 = function(tvec, parameter_vector) {
      stop("CGF$K2 is not initialized.", call. = FALSE)
    },
    K3operator = function(tvec, parameter_vector, v1, v2, v3) {
      stop("CGF$K3operator is not initialized.", call. = FALSE)
    },
    K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) {
      stop("CGF$K4operator is not initialized.", call. = FALSE)
    },

    initialize = function(
      K, K1, K2, K3operator, K4operator,
      analytic_tvec_hat = NULL,
      rsim = NULL,
      op_name = "UnnamedOperation",
      ...
    ) {
      # Make supplied functions behave like R6 methods (access to self/private).
      # Keeps the original scope via parent.env().
      as_method <- function(m) {
        if (!is.function(m)) return(m)
        env_with_self <- new.env(parent = environment(fun = m), size = 2L, hash = FALSE)
        assign("self", value = self, envir = env_with_self)
        assign("private", value = private, envir = env_with_self)
        environment(m) <- env_with_self
        m
      }

      # Wrap a user-supplied function to match the formal arguments of a target
      # method, while still calling the user function with its preferred argument
      # names (mapped positionally) to avoid unused-argument errors.
      wrap_like <- function(target_fun, user_fun) {
        if (!is.function(target_fun)) {
          stop("Internal error: target method is not a function.", call. = FALSE)
        }
        if (!is.function(user_fun)) {
          stop("Internal error: supplied override is not a function.", call. = FALSE)
        }

        target_formals <- formals(target_fun)
        target_names <- names(target_formals)

        user_names <- names(formals(user_fun))
        if (is.null(user_names)) user_names <- character(0)

        dots_pos <- match("...", user_names, nomatch = length(user_names) + 1L)
        n_map <- min(length(target_names), dots_pos - 1L, length(user_names))

        call_names <- target_names
        if (n_map > 0L) call_names[seq_len(n_map)] <- user_names[seq_len(n_map)]

        call_args <- lapply(target_names, as.name)
        names(call_args) <- call_names
        call_expr <- as.call(c(list(quote(user_fun)), call_args))

        env <- new.env(parent = environment(user_fun), size = 1L, hash = FALSE)
        env$user_fun <- user_fun
        eval(call("function", as.pairlist(target_formals), call_expr), env)
      }

      if (!is.function(K)) stop("'K' must be a function.", call. = FALSE)
      if (!is.function(K1)) stop("'K1' must be a function.", call. = FALSE)
      if (!is.function(K2)) stop("'K2' must be a function.", call. = FALSE)
      if (!is.function(K3operator)) stop("'K3operator' must be a function.", call. = FALSE)
      if (!is.function(K4operator)) stop("'K4operator' must be a function.", call. = FALSE)

      required <- list(K = K, K1 = K1, K2 = K2, K3operator = K3operator, K4operator = K4operator)
      for (nm in names(required)) {
        user_fun <- as_method(required[[nm]])
        target_fun <- self[[nm]]
        wrapped <- wrap_like(target_fun, user_fun)
        unlockBinding(nm, self)
        self[[nm]] <- wrapped
        lockBinding(nm, self)
      }

      if (!is.null(analytic_tvec_hat)) {
        if (!is.function(analytic_tvec_hat)) stop("'analytic_tvec_hat' must be NULL or a function.", call. = FALSE)
        self$has_analytic_tvec_hat <- TRUE
        private$analytic_tvec_hat_func <- as_method(analytic_tvec_hat)
        self$analytic_tvec_hat <- function(x, parameter_vector) {
          if (!is.numeric(x)) stop("'x' must be numeric.", call. = FALSE)
          if (any(!is.finite(x))) stop("'x' must be finite.", call. = FALSE)
          private$analytic_tvec_hat_func(x, parameter_vector)
        }
      } else {
        self$has_analytic_tvec_hat <- FALSE
        private$analytic_tvec_hat_func <- NULL
        self$analytic_tvec_hat <- NULL
      }

      if (!is.null(rsim)) {
        if (!is.function(rsim)) stop("'rsim' must be NULL or a function.", call. = FALSE)
        self$has_rsim <- TRUE
        private$rsim_func <- as_method(rsim)
      } else {
        self$has_rsim <- FALSE
        private$rsim_func <- NULL
      }

      if (!is.character(op_name)) stop("'op_name' must be of type character", call. = FALSE)
      self$call_history <- if (!is.null(self$call_history)) c(self$call_history, op_name) else op_name

      # Extra named method overrides in ...
      extra_args <- list(...)
      n_extra <- length(extra_args)
      extra_names <- names(extra_args)
      if (n_extra > 0L && is.null(extra_names)) {
        extra_names <- rep("", n_extra)
      }

      keep_non_null <- !vapply(extra_args, is.null, logical(1))
      if (any(!keep_non_null)) {
        extra_args <- extra_args[keep_non_null]
        extra_names <- extra_names[keep_non_null]
      }
      if (length(extra_args) > 0L && any(!nzchar(extra_names))) {
        warning("Unnamed entries in '...' are ignored. Please provide named overrides.", call. = FALSE)
      }

      override_name <- function(n, e) {
        m <- extra_args[[n]]
        if (!is.function(m)) {
          stop("Override for '", n, "' must be a function (or NULL to keep default).", call. = FALSE)
        }
        m <- as_method(m)
        target_fun <- get(n, envir = e, inherits = FALSE)
        if (!is.function(target_fun)) {
          stop("Cannot override non-function member '", n, "'.", call. = FALSE)
        }
        m <- wrap_like(target_fun, m)
        unlockBinding(n, e)
        assign(n, m, envir = e)
        lockBinding(n, e)
      }

      private_subset <- (nzchar(extra_names)) & (extra_names %in% names(CGF_private_defaults))
      if (any(private_subset)) {
        lapply(extra_names[private_subset], override_name, e = private)
      }

      public_subset <- (nzchar(extra_names)) & (extra_names %in% names(CGF_public_defaults))
      if (any(public_subset)) {
        lapply(extra_names[public_subset], override_name, e = self)
      }

      additional_subset <- (nzchar(extra_names)) &
        !(extra_names %in% names(CGF_private_defaults)) &
        !(extra_names %in% names(CGF_public_defaults))
      if (any(additional_subset)) {
        extras <- lapply(extra_args[additional_subset], as_method)
        self$additional_methods <- modifyList(self$additional_methods, extras)
      }
    },

    print = function(...) {
      cat("<CGF Object>\n")
      if (!is.null(self$call_history)) {
        mapping_str <- paste(self$call_history, collapse = " -> ")
        cat("Used:", mapping_str, "\n")
      }
      cat("Class hierarchy:", paste(class(self), collapse = " -> "), "\n")
      invisible(self)
    },

		    compute.spa.negll = function(parameter_vector,
		                                 observed.data,
		                                 tvec.hat = NULL,
		                                 gradient = FALSE,
	                                 hessian  = FALSE,
	                                 spa_method = "standard",
	                                 ...) {
	      compute.spa.negll(
	        cgf              = self,
	        parameter_vector = parameter_vector,
	        observed.data    = observed.data,
	        tvec.hat         = tvec.hat,
	        gradient         = gradient,
	        hessian          = hessian,
	        spa_method       = spa_method,
	        ...
	      )
	    }
  ))
)


# ------------------------------------------------------------------------
# FACTORY FUNCTION: createCGF
# ------------------------------------------------------------------------
#' Create a CGF object from user-defined functions
#'
#' @description
#' This creates an object of type `CGF` using user-supplied functions. You supply
#' the five essential methods (`K`, `K1`, `K2`, `K3operator`, `K4operator`) plus
#' any optional overrides (e.g., `tilting_exponent` or `neg_ll`), and it returns
#' a `CGF` instance.
#'
#' @param K A function `K(tvec, parameter_vector) -> numeric scalar`.
#' @param K1 A function `K1(tvec, parameter_vector) -> numeric vector`.
#' @param K2 A function `K2(tvec, parameter_vector) -> numeric matrix`.
#' @param K3operator A function implementing the third-order operator.
#' @param K4operator A function implementing the fourth-order operator.
#'
#' @param ineq_constraint Optional function for inequality constraints.
#' @param analytic_tvec_hat Optional function for an analytic solution
#'   of the saddlepoint equation. If provided, call it via `cgf$analytic_tvec_hat(x, param)`.
#' @param rsim Optional simulation method. A function of the form
#'   `function(n, vector_length, parameter_vector, tvec = NULL, ...)` returning
#'   a numeric vector of length `n * vector_length` or a `vector_length x n` matrix.
#'   If supplied, the resulting CGF exposes `$rsim()` and `$has_rsim`.
#' @param op_name A descriptive label for the CGF object/operation. Default is "UnnamedOperation".
#'
#' @param tilting_exponent (optional) Overriding function for the tilting exponent.
#' @param neg_ll (optional) Overriding function for the negative log-likelihood.
#' @param func_T (optional) Overriding function for the first-order correction term.
#' @param K2_solve,logdetK2 (optional) Overriding numerical helper methods.
#' @param K2operator,K2operatorAK2AT,K4operatorAABB,K3K3operatorAABBCC,K3K3operatorABCABC (optional) Overriding operator methods.
#' @param K4operatorAABB_factored,K3K3operatorAABBCC_factored,K3K3operatorABCABC_factored (optional) Overriding factored-operator methods.
#' @param ... Additional named methods or overrides.
#'
#' @return An object of class `CGF`.
#' @export
createCGF <- function(K, K1, K2, K3operator, K4operator,
                      ineq_constraint = NULL,
                      analytic_tvec_hat = NULL,
                      rsim = NULL,
                      op_name = "UnnamedOperation",
                      tilting_exponent = NULL,
                      neg_ll = NULL,
                      func_T = NULL,
                      K2_solve = NULL,
                      logdetK2 = NULL,
                      K4operatorAABB = NULL,
                      K3K3operatorAABBCC = NULL,
                      K3K3operatorABCABC = NULL,
                      K4operatorAABB_factored = NULL,
                      K3K3operatorAABBCC_factored = NULL,
                      K3K3operatorABCABC_factored = NULL,
                      K2operator = NULL,
                      K2operatorAK2AT = NULL,
                      ...
) {
  user_optional_methods <- list(
    ineq_constraint             = ineq_constraint,
    tilting_exponent            = tilting_exponent,
    neg_ll                      = neg_ll,
    func_T                      = func_T,
    K2_solve                    = K2_solve,
    logdetK2                    = logdetK2,
    K4operatorAABB              = K4operatorAABB,
    K3K3operatorAABBCC          = K3K3operatorAABBCC,
    K3K3operatorABCABC          = K3K3operatorABCABC,
    K4operatorAABB_factored     = K4operatorAABB_factored,
    K3K3operatorAABBCC_factored = K3K3operatorAABBCC_factored,
    K3K3operatorABCABC_factored = K3K3operatorABCABC_factored,
    K2operator                  = K2operator,
    K2operatorAK2AT             = K2operatorAK2AT
  )

  additional_methods <- list(...)
  if (length(additional_methods) > 0L && is.null(names(additional_methods))) {
    warning("Unnamed entries in '...' are ignored. Please provide named overrides.", call. = FALSE)
    additional_methods <- list()
  }

  all_optional_methods <- modifyList(user_optional_methods, additional_methods)
  all_optional_methods <- all_optional_methods[!vapply(all_optional_methods, is.null, logical(1))]

  do.call(CGF$new, c(
    list(
      K                 = K,
      K1                = K1,
      K2                = K2,
      K3operator        = K3operator,
      K4operator        = K4operator,
      analytic_tvec_hat = analytic_tvec_hat,
      rsim              = rsim,
      op_name           = op_name
    ),
    all_optional_methods
  ))
}
