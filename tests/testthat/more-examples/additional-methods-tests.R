# Sanity checks for `additional_methods`.
#
# This is a standalone script (not run by testthat) intended for manual use.
# Run from the package root with:
#   R -q -f tests/testthat/more-examples/additional-methods-tests.R

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(".")
} else {
  library(saddlepoint)
}

assert_equal <- function(x, y, tol = 1e-12) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tol)))
}


# -------------------------------------------------------------------------
# Example 1: functions get `self` / `private` injected
# -------------------------------------------------------------------------

toy <- createCGF(
  K = function(tvec, parameter_vector) {
    mu <- parameter_vector[1]
    sigma <- parameter_vector[2]
    sum(mu * tvec + 0.5 * sigma^2 * tvec^2)
  },
  K1 = function(tvec, parameter_vector) {
    mu <- parameter_vector[1]
    sigma <- parameter_vector[2]
    mu + sigma^2 * tvec
  },
  K2 = function(tvec, parameter_vector) {
    sigma <- parameter_vector[2]
    diag(rep(sigma^2, length(tvec)))
  },
  K3operator = function(tvec, parameter_vector, v1, v2, v3) 0,
  K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) 0,
  op_name = "ToyNormal",
  te_via_private = function(tvec, parameter_vector) private$tilting_exponent(tvec, parameter_vector),
  k_via_self = function(tvec, parameter_vector) self$K(tvec, parameter_vector),
  tag = "demo"
)

theta <- c(mu = 1.25, sigma = 0.7)
tvec <- c(-0.2, 0.1, 0.0, 0.3)

stopifnot(is.function(toy$additional_methods$k_via_self))
stopifnot(is.function(toy$additional_methods$te_via_private))
stopifnot(identical(toy$additional_methods$tag, "demo"))

assert_equal(
  as.numeric(toy$additional_methods$k_via_self(tvec, theta)),
  as.numeric(toy$K(tvec, theta))
)

te_ref <- toy$.get_private_method("tilting_exponent")(tvec, theta)
assert_equal(
  as.numeric(toy$additional_methods$te_via_private(tvec, theta)),
  as.numeric(te_ref)
)


# -------------------------------------------------------------------------
# Example 2: lexical scope is preserved
# -------------------------------------------------------------------------

# `multiplier` lives in the scope where `scaled_K` is defined. After we
# "method-ize" additional_methods (injecting `self` / `private` via an enclosing
# environment), this value should still be visible via normal lexical scoping.
multiplier <- 3.0

toy2 <- createCGF(
  K = function(tvec, parameter_vector) sum(parameter_vector[1] * tvec),
  K1 = function(tvec, parameter_vector) rep(parameter_vector[1], length(tvec)),
  K2 = function(tvec, parameter_vector) diag(0, nrow = length(tvec)),
  K3operator = function(tvec, parameter_vector, v1, v2, v3) 0,
  K4operator = function(tvec, parameter_vector, v1, v2, v3, v4) 0,
  op_name = "ToyLinear",
  scaled_K = function(tvec, parameter_vector) multiplier * self$K(tvec, parameter_vector)
)

theta2 <- c(2.0)
tvec2 <- c(0.1, -0.2, 0.3)

assert_equal(
  as.numeric(toy2$additional_methods$scaled_K(tvec2, theta2)),
  multiplier * as.numeric(toy2$K(tvec2, theta2))
)


# -------------------------------------------------------------------------
# Example 3: works when passed via sumOfIndependentCGF(...)
# -------------------------------------------------------------------------

lambda1 <- 2
lambda2 <- 5

cg1 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda1), iidReps = 1)
cg2 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda2), iidReps = 1)
cgf_list <- list(cg1, cg2)

cg_sum <- sumOfIndependentCGF(
  cgf_list,
  K_terms = function(tvec, parameter_vector) {
    vapply(cgf_list, function(cg) cg$K(tvec, parameter_vector), numeric(1))
  },
  K_via_self = function(tvec, parameter_vector) self$K(tvec, parameter_vector)
)

theta_dummy <- 0
t0 <- 0.2

assert_equal(
  as.numeric(sum(cg_sum$additional_methods$K_terms(t0, theta_dummy))),
  as.numeric(cg_sum$K(t0, theta_dummy))
)

assert_equal(
  as.numeric(cg_sum$additional_methods$K_via_self(t0, theta_dummy)),
  as.numeric(cg_sum$K(t0, theta_dummy))
)

cat("OK: additional_methods sanity checks passed.\n")
