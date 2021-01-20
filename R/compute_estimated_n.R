compute_estimated_n_raw <- function(est, theo, only_diag = FALSE){
  if(only_diag){
    x <- diag(theo)
    y <- diag(est)
  } else {
    x <- triangle2vector(theo, diag = TRUE)
    y <- triangle2vector(est, diag = TRUE)
  }

  return(lm(x ~ 0 + y)$coef)
}


compute_estimated_n <- function(dt, only_diag = TRUE){
  dt <- convert_corr_array_to_data_matrix_test(dt)
  est <- t(dt) %*% dt / nrow(dt)
  theo <- corrmat_covariance_from_dt(dt)
  out <- compute_estimated_n_raw(est = est, theo = theo, only_diag = only_diag)
  return(out)
}
