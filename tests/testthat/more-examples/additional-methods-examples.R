# Examples: storing extra helpers on a CGF via `additional_methods`
#
# Anything passed to `createCGF(..., ...)` that is *not* a known CGF method name
# (public or private) will be stored in `cgf$additional_methods`.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(".")
} else {
  library(saddlepoint)
}

# Example 1: attach helpers to a composite CGF (sum of independent summands)
lambda1 <- 2
lambda2 <- 5

cg1 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda1), iidReps = 1)
cg2 <- PoissonModelCGF(lambda = adaptor(fixed_param = lambda2), iidReps = 1)
cg_list <- list(cg1, cg2)

cg_sum <- sumOfIndependentCGF(
  cg_list,
  # Helper: inspect per-summand contributions (uses lexical scoping: cg_list)
  K_by_component = function(tvec, parameter_vector) {
    vapply(cg_list, function(cg) cg$K(tvec, parameter_vector), numeric(1))
  },

  # Helper: anything in the calling scope can be captured too (here: cg_list)
  mean_at0 = function(parameter_vector) {
    sum(vapply(cg_list, function(cg) cg$K1(tvec = 0, parameter_vector = parameter_vector), numeric(1)))
  },

  # Non-function extras can be stored too
  component_names = c("Y1", "Y2")
)

theta_dummy <- 0
t0 <- 0.2

cg_sum$K1(tvec = 0, parameter_vector = theta_dummy) # 7 (= lambda1 + lambda2)
cg_sum$additional_methods$K_by_component(tvec = t0, parameter_vector = theta_dummy)
cg_sum$additional_methods$mean_at0(parameter_vector = theta_dummy)
cg_sum$additional_methods$component_names


# Example 2: you can also add helpers after construction (no `self$...` needed)
cg_sum$additional_methods$K_at0 <- function(parameter_vector) cg_sum$K(tvec = 0, parameter_vector = parameter_vector)
cg_sum$additional_methods$K_at0(theta_dummy)
