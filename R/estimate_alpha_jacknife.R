# todo: beautify
estimate_alpha_jacknife <- function(
  healthy_dt, sick_dt, dim_alpha = 1,
  alpha0 = NULL, theta0 = NULL,
  linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(iter_config = list(min_loop = 0)), cov_config = list(),
  return_gee = FALSE, jack_healthy = TRUE,
  bias_correction = FALSE, early_stop = FALSE,
  verbose = TRUE, ncores = 1
){

  mclapply_ <- if(verbose) pbmcapply::pbmclapply else parallel::mclapply

  apply_fun <- function(i, boot_dt){
    if(boot_dt == 'sick'){
      sick_dt_ <- sick_dt[-i,]
      healthy_dt_ <- healthy_dt
    } else if(boot_dt == 'healthy'){
      sick_dt_ <- sick_dt
      healthy_dt_ <- healthy_dt[-i,]
    }

    weight_matrix <- corrmat_covariance_from_dt(sick_dt_)

    cov_model <- inner_optim_loop(
      healthy_dt = healthy_dt_, sick_dt = sick_dt_,
      alpha0 = alpha0, theta0 = theta0, dim_alpha = dim_alpha,
      weight_matrix = weight_matrix, linkFun = linkFun,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
      early_stop = early_stop, verbose = FALSE
    )

    if(bias_correction) cov_model$alpha <- cov_model$alpha - median(cov_model$alpha) + linkFun$null_value

    gee_out <- if(return_gee){
      triangle2vector(
        compute_gee_variance(
          cov_obj = cov_model,
          healthy_dt = healthy_dt_,
          sick_dt = sick_dt_,
          est_mu = TRUE
        ),
        diag = TRUE
      )
    } else NA

    return(list(
      theta = cov_model$theta,
      alpha = cov_model$alpha,
      convergence = tail(cov_model$convergence, 1),
      gee_var = gee_out
    ))
  }

  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }

  healthy_dt <- convert_corr_array_to_data_matrix(healthy_dt)
  sick_dt <- convert_corr_array_to_data_matrix(sick_dt)

  if(is.null(alpha0) | is.null(theta0)){
    ini_model <- estimate_alpha(
      healthy_dt = healthy_dt,
      sick_dt = sick_dt,
      dim_alpha = dim_alpha,
      linkFun = linkFun,
      model_reg_config = model_reg_config,
      matrix_reg_config = matrix_reg_config,
      raw_start = TRUE,
      iid_config = iid_config,
      cov_config = cov_config,
      bias_correction = bias_correction,
      early_stop = early_stop,
      verbose = FALSE
    )
    index <- if(length(ini_model$steps) > 3) (length(ini_model$steps) - 3) else (length(ini_model$steps) - 1)
    ini_model <- ini_model$steps[[index]]
    if(is.null(alpha0)) alpha0 <- as.vector(ini_model$alpha)
    if(is.null(theta0)) theta0 <- ini_model$theta
  }

  if(verbose) cat('\njacknifing Sick Observations...\n')
  cov_obj_sick <- mclapply_(seq_len(nrow(sick_dt)), apply_fun, boot_dt = 'sick', mc.cores = ncores)
  cov_obj_sick_t <- purrr::transpose(cov_obj_sick)

  theta <- do.call(rbind, cov_obj_sick_t$theta)
  alpha <- do.call(rbind, lapply(cov_obj_sick_t$alpha, as.vector))
  convergence <- do.call(c, cov_obj_sick_t$convergence)
  gee_var <- do.call(rbind, cov_obj_sick_t$gee_var)

  if(jack_healthy){
    if(verbose) cat('\njacknifing Healthy Observations...\n')
    cov_obj_healthy <- mclapply_(seq_len(nrow(healthy_dt)), apply_fun, boot_dt = 'healthy', mc.cores = ncores)
    cov_obj_healthy_t <- purrr::transpose(cov_obj_healthy)

    theta_h <- do.call(rbind, cov_obj_healthy_t$theta)
    alpha_h <- do.call(rbind, lapply(cov_obj_healthy_t$alpha, as.vector))
    convergence_h <- do.call(c, cov_obj_healthy_t$convergence)
    gee_var_h <- do.call(rbind, cov_obj_healthy_t$gee_var)

    theta <- rbind(theta_h, theta)
    alpha <- rbind(alpha_h, alpha)
    convergence <- c(convergence_h, convergence)
    gee_var <- rbind(gee_var, gee_var_h)
    is_sick <- c(rep(0, nrow(healthy_dt)), rep(1, nrow(sick_dt)))
  }

  return(list(
    theta = theta,
    alpha = alpha,
    gee_var = gee_var,
    convergence = convergence,
    is_sick = is_sick,
    linkFun = linkFun
  ))
}
