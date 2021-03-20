default_model_reg_config <- list(lambda = 0, lp = 2)
default_matrix_reg_config <- list(method = 'constant', const = 0)
default_iterations_config <- list(maxit = 50, reltol = 1e-06, minit = 3)
default_optim_config <- list(method = "BFGS", reltol = 1e-06, log_optim = FALSE)


configs <- list(
  model_reg = default_model_reg_config,
  matrix_reg = default_matrix_reg_config,
  iterations = default_iterations_config,
  optim = default_optim_config
)


#' Default Configurations for Model Estimation
#'
#' A list of default configurations for \link[corrpops]{estimate_model}. running the function will print the defaults.
#' - model_reg: a list with params lambda, lp for regularization of alpha. the loss function has an element of reg_lambda * sum(abs(alpha - LinkFunc$null_value)^reg_p) ^ (1/reg_p)
#' - matrix_reg: parameters passed on to \link[corrpops]{regularize_matrix}. the weighting matrix (covariance matrix of correlations) is regularized accordingly
#' - iterations: parameters used for the iterations' stopping rule. the definitions of 'maxit', 'reltol', 'abstol' are the same is in \link[stats]{optim}. iterations will stop only if there were 'minit' optim iterations with convergence=0 in a row.
#' - optim: possible parameters to be passed to \link[stats]{optim}, namely 'method', 'reltol', 'abstol'. also, log_optim is a boolean - whether to save the result of the call to \link[stats]{optim} in each iteration.
#' @param index which configurations to print. if missing, print all. can be one or more of the aforementioned.
#' @export
configurations <- function(index)
  if(missing(index)){
    print(configs)
  } else {
    print(configs[index])
  }
