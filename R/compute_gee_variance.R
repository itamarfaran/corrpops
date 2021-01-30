#' Compute the Jacobian by Alpha parameter
#'
#' Helper function to calulcate the Jacobian of the expected value by Alpha:
#' \eqn{J_{\alpha}\left(g\left(\Theta,\alpha\right)\right)_{ij}=
#' \frac{\partial g_{i}\left(\Theta,\alpha\right)}{\partial\alpha_{j}}}
#'
#' @param group can be either 'diagnosed' or 'control'
#' @param alpha the vector to calculate the Jacobian on
#' @param control_dt array of correlation matrices of control group, vectorized and organized in a data matrix
#' @param diagnosed_dt same as control_dt but of the diagnosed group
#' @param d the number of columns in alpha, default 1
#' @param linkFunc the link funtion $g$ to use
#' @return the Jacobian matrix
#'
compute_mu_alpha_jacobian <- function(group, alpha, control_dt, diagnosed_dt, d = 1, LinkFunc)
{
   if(group == 'diagnosed'){
     func <- function(A){
       out <- triangle2vector(LinkFunc$func(
         t = theta_of_alpha(A, control_dt, diagnosed_dt, LinkFunc = LinkFunc, d = d),
         a = A,
         d = d
       ))
       return(out)
     }
   } else if(group == 'control'){
     func <- function(A){
       out <- theta_of_alpha(A, control_dt, diagnosed_dt, LinkFunc = LinkFunc, d = d)
       return(out)
     }
   }
  out <- jacobian(func = func, x = alpha)
  return(out)
}


#' Compute Estimates' Variance Using a GEE Framework
#'
#' Estimate the estimators' variance using a Generelized Estimation Equations framwork
#'
#' @param mod an object created by \link[corrfuncs]{estimate_model}
#' @param control_dt array of correlation matrices of control group. can be passed as a matrix, each row is a vectorized correlation matrix
#' @param diagnosed_dt same as control_dt but of the diagnosed group
#' @param est_mu whether to use the expected value estimated in the model or the sample average. default TRUE.
#' @param reg_lambda todo: drop, should inherit from mod
#' @param reg_p todo: drop, should inherit from mod
#' @return fisher information matrix of alpha
#'
#' @export
#'
compute_gee_variance <- function(mod,
                                 control_dt, diagnosed_dt,
                                 est_mu = TRUE, reg_lambda = 0, reg_p = 2)
{
  inner <- function(group){
    p <- 0.5 + sqrt(1 + 8*ncol(data))/2
    d <- length(mod$alpha)/p

    control_data <- convert_corr_array_to_data_matrix(control_dt)
    diagnosed_data <- convert_corr_array_to_data_matrix(diagnosed_dt)
    data <- if(group == 'control') control_data else diagnosed_data

    jacobian <- compute_mu_alpha_jacobian(
      group = group,
      alpha = mod$alpha,
      control_dt = control_data,
      diagnosed_dt = diagnosed_data,
      d = d,
      LinkFunc = mod$LinkFunc)

    if(est_mu){
      if(group == 'control'){
        expected_value <- mod$theta
      } else {
        expected_value <- triangle2vector(mod$LinkFunc$func(
          t = mod$theta,
          a = mod$alpha,
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

  control_lst <- inner('control')
  diagnosed_lst <- inner('diagnosed')

  I0 <- compute_gee_raw('I0', control_lst) + compute_gee_raw('I0', diagnosed_lst)

  I1 <- compute_gee_raw('I1', control_lst) + compute_gee_raw('I1', diagnosed_lst)
  solve_I0 <- solve(I0)

  out <- solve_I0 %*% I1 %*% solve_I0
  return(out)
}
