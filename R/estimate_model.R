estimate_model <- function(
  control_dt, diagnosed_dt, dim_alpha = 1,
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

  control_dt <- convert_corr_array_to_data_matrix(control_dt)
  diagnosed_dt <- convert_corr_array_to_data_matrix(diagnosed_dt)

  alpha0 <- theta0 <- NULL
  if(!raw_start){
    iid_model <- optimiser(
      control_dt = control_dt, diagnosed_dt = diagnosed_dt,
      alpha0 = NULL, theta0 = NULL,
      weight_matrix = NULL, dim_alpha = dim_alpha, LinkFunc = LinkFunc,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = iid_config$iter_config, optim_config = iid_config$optim_config,
      early_stop = early_stop, verbose = verbose[1]
    )
    alpha0 <- iid_model$alpha
    theta0 <- iid_model$theta
  }

  weight_matrix <- corrmat_covariance_from_dt(diagnosed_dt)

  cov_model <- optimiser(
    control_dt = control_dt, diagnosed_dt = diagnosed_dt,
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
