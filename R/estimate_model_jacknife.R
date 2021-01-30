estimate_model_jacknife <- function(
  control_dt, diagnosed_dt, dim_alpha = 1, alpha0 = NULL, theta0 = NULL,
  LinkFunc = LinkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(iter_config = list(min_loop = 0)), cov_config = list(),
  return_gee = FALSE, jack_control = TRUE,
  bias_correction = FALSE, early_stop = FALSE,
  verbose = FALSE, ncores = 1
){
  if(ncores > 1){
    if(verbose){
      lapply_ <- pbmcapply::pbmclapply
    } else {
      lapply_ <- parallel::mclapply
    }
  } else {
     if(verbose){
      lapply_ <- pbapply::pblapply
    } else {
      lapply_ <- lapply
    }
  }

  apply_func <- function(i, boot_dt){
    if(boot_dt == 'diagnosed'){
      diagnosed_dt_ <- diagnosed_dt[-i,]
      control_dt_ <- control_dt
    } else if(boot_dt == 'control'){
      diagnosed_dt_ <- diagnosed_dt
      control_dt_ <- control_dt[-i,]
    }

    weight_matrix <- corrmat_covariance_from_dt(diagnosed_dt_)
    cov_model <- optimiser(
      control_dt = control_dt_, diagnosed_dt = diagnosed_dt_,
      alpha0 = alpha0, theta0 = theta0, dim_alpha = dim_alpha,
      weight_matrix = weight_matrix, LinkFunc = LinkFunc,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
      early_stop = early_stop, verbose = FALSE
    )

    if(bias_correction)
      cov_model$alpha <- cov_model$alpha - median(cov_model$alpha) + LinkFunc$null_value

    gee_out <- NA
    if(return_gee){
      gee_out <- triangle2vector(compute_gee_variance(
        mod = cov_model,
        control_dt = control_dt_,
        diagnosed_dt = diagnosed_dt_,
        est_mu = TRUE),
        diag = TRUE)
    }

    out <- list(
      theta = cov_model$theta,
      alpha = cov_model$alpha,
      convergence = tail(cov_model$convergence, 1),
      gee_var = gee_out
    )

    return(out)
  }

  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }

  control_dt <- convert_corr_array_to_data_matrix(control_dt)
  diagnosed_dt <- convert_corr_array_to_data_matrix(diagnosed_dt)

  if(is.null(alpha0) | is.null(theta0)){
    ini_model <- estimate_alpha(
      control_dt = control_dt,
      diagnosed_dt = diagnosed_dt,
      dim_alpha = dim_alpha,
      LinkFunc = LinkFunc,
      model_reg_config = model_reg_config,
      matrix_reg_config = matrix_reg_config,
      raw_start = TRUE,
      iid_config = iid_config,
      cov_config = cov_config,
      bias_correction = bias_correction,
      early_stop = early_stop,
      verbose = FALSE
    )

    ini_model_steps <- length(ini_model$steps)
    index <- ini_model_steps - 1
    if(ini_model_steps > 3)
      index <- ini_model_steps - 3

    ini_model <- ini_model$steps[[index]]
    if(is.null(alpha0))
      alpha0 <- as.vector(ini_model$alpha)
    if(is.null(theta0))
      theta0 <- ini_model$theta
  }

  if(verbose)
    cat('\njacknifing diagnosed observations...\n')

  if(ncores > 1){
    mod_diagnosed <- lapply_(seq_len(nrow(diagnosed_dt)), apply_func, boot_dt = 'diagnosed', mc.cores = ncores)
  } else {
    mod_diagnosed <- lapply_(seq_len(nrow(diagnosed_dt)), apply_func, boot_dt = 'diagnosed')
  }

  mod_diagnosed_t <- purrr::transpose(mod_diagnosed)

  theta <- do.call(rbind, mod_diagnosed_t$theta)
  alpha <- do.call(rbind, lapply(mod_diagnosed_t$alpha, as.vector))
  convergence <- do.call(c, mod_diagnosed_t$convergence)
  gee_var <- do.call(rbind, mod_diagnosed_t$gee_var)

  if(jack_control){
    if(verbose)
      cat('\njacknifing control observations...\n')
    if(ncores > 1){
      mod_control <- lapply_(seq_len(nrow(control_dt)), apply_func, boot_dt = 'control', mc.cores = ncores)
    } else {
      mod_control <- lapply_(seq_len(nrow(control_dt)), apply_func, boot_dt = 'control')
    }

    mod_control_t <- purrr::transpose(mod_control)

    theta_h <- do.call(rbind, mod_control_t$theta)
    alpha_h <- do.call(rbind, lapply(mod_control_t$alpha, as.vector))
    convergence_h <- do.call(c, mod_control_t$convergence)
    gee_var_h <- do.call(rbind, mod_control_t$gee_var)

    theta <- rbind(theta_h, theta)
    alpha <- rbind(alpha_h, alpha)
    convergence <- c(convergence_h, convergence)
    gee_var <- rbind(gee_var, gee_var_h)
    is_diagnosed <- c(rep(0, nrow(control_dt)), rep(1, nrow(diagnosed_dt)))
  }

  output <- list(
    theta = theta,
    alpha = alpha,
    gee_var = gee_var,
    convergence = convergence,
    is_diagnosed = is_diagnosed,
    LinkFunc = LinkFunc
  )

  return(output)
}
