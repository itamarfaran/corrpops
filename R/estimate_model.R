#' Estimate Model
#'
#' todo: explain
#'
#' @param control_arr array of control group correlation matrices. either an array or data matrix form
#' @param diagnosed_arr array of diagnosed group correlation matrices. either an array or data matrix form
#' @param dim_alpha the number of columns in alpha. default 1
#' @param LinkFunc a list of function. must include func, inverse, rev_func and null_value. @seealso LinkFunc
#' @param model_reg_config list of configurations for the model regularization seealso
#' @param matrix_reg_config list of configurations for the covariance matrix's regularization seealso
#' @param iid_config list of configurations for the optimization of the model with identity matrix covariance matrix seealso
#' @param cov_config list of configurations for the optimization of the model with specified covariance matrix seealso
#' @param raw_start if true, don't optimize with the identity matrix before optimizing with a specified covariance matrix
#' @param bias_correction if true, correct the estimates to the median: a' = a - med(a) + null_value
#' @param early_stop if true, stop the optimization of the joint loss function (of theta and alpha) didn't decrease.
#' @param verbose if true, print status to console
#' @return a list of todo
#'
#' @export
#'
estimate_model <- function(
  control_arr, diagnosed_arr, dim_alpha = 1,
  LinkFunc = LinkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(), cov_config = list(),
  raw_start = FALSE, bias_correction = FALSE, early_stop = FALSE,
  verbose = TRUE)
{
  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }
  if(length(verbose) == 1) verbose <- rep(verbose, 2)

  control_datamatrix <- convert_corr_array_to_data_matrix(control_arr)
  diagnosed_datamatrix <- convert_corr_array_to_data_matrix(diagnosed_arr)

  alpha0 <- theta0 <- NULL
  if(!raw_start){
    iid_model <- optimiser(
      control_datamatrix = control_datamatrix, diagnosed_datamatrix = diagnosed_datamatrix,
      alpha0 = NULL, theta0 = NULL,
      weight_matrix = NULL, dim_alpha = dim_alpha, LinkFunc = LinkFunc,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = iid_config$iter_config, optim_config = iid_config$optim_config,
      early_stop = early_stop, verbose = verbose[1]
    )
    alpha0 <- iid_model$alpha
    theta0 <- iid_model$theta
  }

  weight_matrix <- corrmat_covariance_from_datamatrix(diagnosed_datamatrix)

  cov_model <- optimiser(
    control_datamatrix = control_datamatrix, diagnosed_datamatrix = diagnosed_datamatrix,
    alpha0 = alpha0, theta0 = theta0,
    weight_matrix = weight_matrix, LinkFunc = LinkFunc,
    model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
    iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
    early_stop = early_stop, verbose = verbose[2]
  )

  if(bias_correction)
    cov_model$alpha <- cov_model$alpha - median(cov_model$alpha) + LinkFunc$null_value

  return(cov_model)
}
