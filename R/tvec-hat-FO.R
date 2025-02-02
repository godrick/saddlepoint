



# Helper function to create dfdu_solve_fn with access to 'cgf'
create_tvec_hat_K2_solve_fn <- function(cgf) {
  function(tvec, theta, w) {
    K2_result <- cgf$K2(tvec, theta)
    solve(K2_result, w)
  }
}


# choose_spa_function <- function(method, cgf) {
#   # Allowed "methods" or "types" of SPA functions
#   # can expand in future => when non-gaussian methods are added
#   # to expand just add them to allowed_methods and create another branch in the if-else
#   allowed_methods <- c(
#     "negll_standard",  # standard spa negative log-likelihood
#     "negll_zeroth",    # zeroth-order negative log-likelihood
#     "correction_standard",   # correction to the first-order/standard spa to the negative log-likelihood
#     "correction_zeroth"      # correction to the zeroth-order spa to the negative log-likelihood
#   )
#   
#   if (!method %in% allowed_methods) {
#     stop(sprintf(
#       "Unknown method '%s'. Allowed methods: %s",
#       method, paste(allowed_methods, collapse = ", ")
#     ))
#   }
#   
#   # Access the CGF private methods just once here:
#   neg_ll_fun  <- cgf$.get_private_method("neg_ll")             # function(t, theta)
#   tilt_fun    <- cgf$.get_private_method("tilting_exponent")   # function(t, theta)
#   funcT_first <- cgf$.get_private_method("func_T")             # (t, theta)
#   
#   # Return whichever function is appropriate:
#   if (method == "negll_standard") {
#     return(neg_ll_fun)
#     
#   } else if (method == "negll_zeroth") {
#     return(function(tvec, parameter_vector) {- tilt_fun(t, parameter_vector)} )
#     
#   } else if (method == "correction_standard") {
#     return(funcT_first)
#     
#   } else {
#     # method == "corr_zeroth" => the function computing -0.5 * log det(K2)
#     return(
#       function(tvec, par_vec) {
#         K2_val <- cgf$K2(tvec, par_vec)
#         -0.5 * determinant(K2_val, logarithm = TRUE)$modulus
#     })
#   }
# }
# 
# 
# 
# get_spa_taped_fun <- function(param_vec,
#                               observed.data,
#                               cgf,
#                               final_tvec,
#                               method) {
#   
#   # method is one of "negll_standard", "negll_zeroth", "corr_standard", "corr_zeroth"
#   chosen_spa_fn <- choose_spa_function(method = method, cgf = cgf)
#   K2_solve_fn   <- create_tvec_hat_K2_solve_fn(cgf)
#   
#   local_fn <- function(par) {
#     tvec_for_ad <- tvec_hat(
#       theta         = par,
#       current_tvec  = final_tvec,
#       observed.data = observed.data,
#       K1_fn         = cgf$K1,
#       K2_solve_fn   = K2_solve_fn
#     )
#     chosen_spa_fn(tvec_for_ad, par)
#   }
#   
#   MakeTape(f = local_fn, x = param_vec)
# }



# get_tape_res <- function(tape_obj, param_vec, gradient = FALSE, hessian = FALSE) {
#   
#   out_list <- list(
#     vals = tape_obj(param_vec),
#     gradient = NULL,
#     hessian = NULL
#   )
#   
#   if (gradient) { out_list$gradient <- as.vector(tape_obj$jacobian(param_vec)) }
#   
#   if (hessian) {
#     jacfun_obj <- tape_obj$jacfun()
#     out_list$hessian <- jacfun_obj$jacobian(param_vec)
#   }
#   out_list
# }






