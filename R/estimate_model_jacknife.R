#' Jackknife Estimation of Model
#'
#' todo: explain
#'
#' @param control_arr array of control group correlation matrices. either an array or data matrix form
#' @param diagnosed_arr array of diagnosed group correlation matrices. either an array or data matrix form
#' @param dim_alpha the number of columns in alpha. default 1
#' @param alpha0 starting point for alpha in the optimization. if null (the default), will use LinkFunc$null_value
#' @param theta0 starting point for alpha in the optimization. if null (the default), will the average matrix of all subjects
#' @param LinkFunc a list of function. must include func, inverse, rev_func and null_value. see \link[corrpops]{LinkFuncSkeleton}
#' @param infer whether to claculate the jackknife framework covariance matrix. default TRUE
#' @param model_reg_config see \link[corrpops]{configurations}. arguments passed will override the defaults.
#' @param matrix_reg_config see \link[corrpops]{configurations}. arguments passed will override the defaults.
#' @param iid_config list of two lists named 'iter_config' and 'optim_config', for the optimization of the model with identity matrix covariance matrix. see \link[corrpops]{configurations}. arguments passed will override the defaults.
#' @param cov_config list of two lists named 'iter_config' and 'optim_config', for the optimization of the model with a specified covariance matrix. see \link[corrpops]{configurations}. arguments passed will override the defaults.
#' @param return_gee if true, calculate the gee estimate of variance in each jackknife
#' @param jack_control if false, don't jackknife control subjects
#' @param bias_correction if true, correct the estimates to the median: a' = a - med(a) + null_value
#' @param early_stop if true, stop the optimization of the joint loss function (of theta and alpha) didn't decrease.
#' @param early_stop if true, stop the optimization of the joint loss function (of theta and alpha) didn't decrease in the last iteration.
#' @param verbose if true, print status to console
#' @param ncores number of cores to use in parallel
#' @return a list of the following:
#' @return - theta: a matrix, the estimates of theta from each jackknife iteration vectorized
#' @return - alpha: a matrix, the estimates of alpha from each jackknife iteration vectorized
#' @return - gee_var: if return_gee=TRUE, an array of the gee estimate of variance from each jackknife iteration
#' @return - convergence: a matrix consisting of vectors with the convergence in each iteration. see \link[stats]{optim}
#' @return - is_diagnosed: a vector with 0,1 indicating of a control or diagnosed subjet was omitted in the iteration
#' @return - LinkFunc: the link function used same as the parameter LinkFunc
#' @return - regularization: the same as model_reg_config parameter
#' @return - vcov: the covariance matrix of the estimates
#' @seealso \link[corrpops]{estimate_model}
#' @export
#'
# todo: estimate_two_pop_model_jk
# todo: add wrapper with jk variance

estimate_model_jacknife <- function(
  control_arr, diagnosed_arr, dim_alpha = 1, alpha0 = NULL, theta0 = NULL,
  LinkFunc = LinkFunctions$multiplicative_identity, infer=TRUE,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(iter_config = list(minit = 0)), cov_config = list(),
  return_gee = FALSE, jack_control = TRUE,
  bias_correction = FALSE, early_stop = FALSE,
  verbose = TRUE, ncores = 1
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

  apply_func <- function(i, boot_datamatrix){
    if(boot_datamatrix == 'diagnosed'){
      diagnosed_datamatrix_ <- diagnosed_datamatrix[-i,]
      control_datamatrix_ <- control_datamatrix
    } else if(boot_datamatrix == 'control'){
      diagnosed_datamatrix_ <- diagnosed_datamatrix
      control_datamatrix_ <- control_datamatrix[-i,]
    }

    weight_matrix <- corrmat_covariance_from_datamatrix(diagnosed_datamatrix_)
    cov_model <- optimiser(
      control_datamatrix = control_datamatrix_, diagnosed_datamatrix = diagnosed_datamatrix_,
      alpha0 = alpha0, theta0 = theta0, dim_alpha = dim_alpha,
      weight_matrix = weight_matrix, LinkFunc = LinkFunc,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
      early_stop = early_stop, verbose = FALSE
    )

    if(bias_correction)
      cov_model$alpha <- cov_model$alpha - stats::median(cov_model$alpha) + LinkFunc$null_value

    gee_out <- NA
    if(return_gee){
      gee_out <- triangle2vector(compute_gee_variance(
        mod = cov_model,
        control_arr = control_datamatrix_,
        diagnosed_arr = diagnosed_datamatrix_,
        est_mu = TRUE),
        diag = TRUE)
    }

    out <- list(
      theta = cov_model$theta,
      alpha = cov_model$alpha,
      convergence = utils::tail(cov_model$convergence, 1),
      gee_var = gee_out
    )

    return(out)
  }

  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }

  control_datamatrix <- convert_corr_array_to_data_matrix(control_arr)
  diagnosed_datamatrix <- convert_corr_array_to_data_matrix(diagnosed_arr)

  if(is.null(alpha0) | is.null(theta0)){
    ini_model <- estimate_model(
      control_arr = control_datamatrix,
      diagnosed_arr = diagnosed_datamatrix,
      dim_alpha = dim_alpha,
      LinkFunc = LinkFunc,
      infer = FALSE,
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
    mod_diagnosed <- lapply_(seq_len(nrow(diagnosed_datamatrix)), apply_func, boot_datamatrix = 'diagnosed', mc.cores = ncores)
  } else {
    mod_diagnosed <- lapply_(seq_len(nrow(diagnosed_datamatrix)), apply_func, boot_datamatrix = 'diagnosed')
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
      mod_control <- lapply_(seq_len(nrow(control_datamatrix)), apply_func, boot_datamatrix = 'control', mc.cores = ncores)
    } else {
      mod_control <- lapply_(seq_len(nrow(control_datamatrix)), apply_func, boot_datamatrix = 'control')
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
    is_diagnosed <- c(rep(0, nrow(control_datamatrix)), rep(1, nrow(diagnosed_datamatrix)))
  }

  output <- list(
    theta = theta,
    alpha = alpha,
    gee_var = gee_var,
    convergence = convergence,
    is_diagnosed = is_diagnosed,
    LinkFunc = LinkFunc,
    regularization = model_reg_config
  )

  output$vcov <- if(infer) infer_jackknife(output) else NULL

  return(output)
}
