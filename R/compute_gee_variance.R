compute_mu_alpha_jacobian <- function(type, alpha, healthy_dt, sick_dt, d = 1, linkFun){
  func <- if(type == 'sick'){
    function(A) triangle2vector(
      linkFun$FUN(
        t = theta_of_alpha(A, healthy_dt, sick_dt, linkFun = linkFun, d = d),
        a = A,
        d = d
      )
    )
  } else if(type == 'healthy') {
    function(A) theta_of_alpha(A, healthy_dt, sick_dt, linkFun = linkFun, d = d)
  }
  return(
    jacobian(func = func, x = alpha)
  )
}


compute_gee_variance <- function(
  cov_obj, healthy_dt, sick_dt, est_mu = TRUE,
  reg_lambda = 0, reg_p = 2){

  create_list_for_raw_gee <- function(type, cov_obj, healthy_dt, sick_dt, est_mu){
    type <- match.arg(type, c('healthy', 'sick'))
    healthy_data <- convert_corr_array_to_data_matrix_test(healthy_dt)
    sick_data <- convert_corr_array_to_data_matrix_test(sick_dt)
    data <- if(type == 'healthy') healthy_data else sick_data

    p <- 0.5 + sqrt(1 + 8*ncol(data))/2
    d <- length(cov_obj$alpha)/p

    jacobian <- compute_mu_alpha_jacobian(
      type = type,
      alpha = cov_obj$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data,
      d = d,
      linkFun = cov_obj$linkFun)
    expected_value <- if(est_mu){
      if(type == 'healthy') cov_obj$theta else triangle2vector(
        cov_obj$linkFun$FUN(
          t = cov_obj$theta,
          a = cov_obj$alpha,
          d = d
        )
      )
    } else colMeans(data)
    solve_Sigma <- solve(corrmat_covariance_from_dt(data, est_n = T))
    df <- nrow(data) - 1
    # efrons_effective_sample_size(
    # n = nrow(data),
    # efrons_rms_sample(data)
    # )

    out <- list(
      data = data,
      jacobian = jacobian,
      expected_value = expected_value,
      solve_Sigma = solve_Sigma,
      df = df
    )
    return(out)
  }

  compute_gee_raw <- function(type, list_){
    if(type == 'I0'){
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% list_$jacobian
    } else if (type == 'I1'){
      residuals <- list_$data - rep(1, nrow(list_$data)) %o% list_$expected_value
      cov_mat <- t(residuals) %*% residuals / list_$df
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% cov_mat %*% list_$solve_Sigma %*% list_$jacobian
    }
    out <- out*nrow(list_$data)
    return(out)
  }

  healthy_list <- create_list_for_raw_gee('healthy', cov_obj, healthy_dt, sick_dt, est_mu)
  sick_list <- create_list_for_raw_gee('sick', cov_obj, healthy_dt, sick_dt, est_mu)

  I0 <- compute_gee_raw('I0', healthy_list) + compute_gee_raw('I0', sick_list)
  solve_I0 <- solve(I0)
  I1 <- compute_gee_raw('I1', healthy_list) + compute_gee_raw('I1', sick_list)
  res <- solve_I0 %*% I1 %*% solve_I0
  return(res)
}
