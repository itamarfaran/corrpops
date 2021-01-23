compute_mu_alpha_jacobian <- function(group, alpha, healthy_dt, sick_dt, d = 1, LinkFunc)
{
   if(group == 'sick'){
     func <- function(A){
       out <- triangle2vector(LinkFunc$func(
         t = theta_of_alpha(A, healthy_dt, sick_dt, LinkFunc = LinkFunc, d = d),
         a = A,
         d = d
       ))
       return(out)
     }
   } else if(group == 'healthy'){
     func <- function(A){
       out <- theta_of_alpha(A, healthy_dt, sick_dt, LinkFunc = LinkFunc, d = d)
       return(out)
     }
   }
  out <- jacobian(func = func, x = alpha)
  return(out)
}


compute_gee_variance <- function(cov_obj,
                                 healthy_dt, sick_dt,
                                 est_mu = TRUE, reg_lambda = 0, reg_p = 2)
{
  inner <- function(group){
    p <- 0.5 + sqrt(1 + 8*ncol(data))/2
    d <- length(cov_obj$alpha)/p

    healthy_data <- convert_corr_array_to_data_matrix(healthy_dt)
    sick_data <- convert_corr_array_to_data_matrix(sick_dt)
    data <- if(group == 'healthy') healthy_data else sick_data

    jacobian <- compute_mu_alpha_jacobian(
      group = group,
      alpha = cov_obj$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data,
      d = d,
      LinkFunc = cov_obj$LinkFunc)

    if(est_mu){
      if(group == 'healthy'){
        expected_value <- cov_obj$theta
      } else {
        expected_value <- triangle2vector(cov_obj$LinkFunc$func(
          t = cov_obj$theta,
          a = cov_obj$alpha,
          d = d))
      }
    } else {
      expected_value <- colMeans(data)
    }

    solve_sigma <- solve(corrmat_covariance_from_dt(data, est_n = T))
    df <- nrow(data) - 1
    # or:
    # efrons_effective_sample_size(
    # n = nrow(data),
    # efrons_rms_sample(data)
    # )

    out <- list(
      data = data,
      jacobian = jacobian,
      expected_value = expected_value,
      solve_sigma = solve_sigma,
      df = df
    )
    return(out)
  }

  compute_gee_raw <- function(group, lst){
    if(group == 'I0'){
      out <- t(lst$jacobian) %*% lst$solve_sigma %*% lst$jacobian
    } else {
      residuals <- lst$data - rep(1, nrow(lst$data)) %o% lst$expected_value
      cov_mat <- t(residuals) %*% residuals / lst$df
      out <- t(lst$jacobian) %*% lst$solve_sigma %*% cov_mat %*% lst$solve_sigma %*% lst$jacobian
    }
    out <- out * nrow(lst$data)
    return(out)
  }

  healthy_lst <- inner('healthy')
  sick_lst <- inner('sick')

  I0 <- compute_gee_raw('I0', healthy_lst) + compute_gee_raw('I0', sick_lst)

  I1 <- compute_gee_raw('I1', healthy_lst) + compute_gee_raw('I1', sick_lst)
  solve_I0 <- solve(I0)

  out <- solve_I0 %*% I1 %*% solve_I0
  return(out)
}
