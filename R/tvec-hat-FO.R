



# Helper function to create dfdu_solve_fn with access to 'cgf'
create_tvec_hat_dfdu_solve_fn <- function(cgf) {
  function(tvec, theta, w) {
    K2_result <- cgf$K2(tvec, theta)
    solve(K2_result, w)
  }
}
