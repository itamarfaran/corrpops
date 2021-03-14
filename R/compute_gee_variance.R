#' Compute the Jacobian of the link function by Alpha parameter
#'
#' Helper function to calculate the Jacobian of the expected value by Alpha:
#' \eqn{J_{\alpha}\left(g\left(\Theta,\alpha\right)\right)_{ij}=
#' \frac{\partial g_{i}\left(\Theta,\alpha\right)}{\partial\alpha_{j}}}
#'
#' @param group can be either 'diagnosed' or 'control'
#' @param alpha the vector to calculate the Jacobian on
#' @param control_datamatrix array of correlation matrices of control group, vectorized and organized in a data matrix
#' @param diagnosed_datamatrix same as control_datamatrix but of the diagnosed group
#' @param d the number of columns in alpha, default 1
#' @param LinkFunc the link function $g$ to use
#' @return the Jacobian matrix
#'
compute_mu_alpha_jacobian <- function(group, alpha, control_datamatrix, diagnosed_datamatrix, d = 1, LinkFunc)
{
   if(group == 'diagnosed'){
     func <- function(A){
       theta <- theta_of_alpha(A, control_datamatrix, diagnosed_datamatrix, LinkFunc = LinkFunc, d = d)

       out <- triangle2vector(
         LinkFunc$func(
           t = theta,
           a = A,
           d = d
           )
         )
       return(out)
     }
   } else if(group == 'control'){
     func <- function(A){
       theta <- theta_of_alpha(A, control_datamatrix, diagnosed_datamatrix, LinkFunc = LinkFunc, d = d)

       out <- theta
       return(out)
     }
   }
  out <- numDeriv::jacobian(func = func, x = alpha)
  return(out)
}


#' Compute Estimates' Variance Using a GEE Framework
#'
#' Estimate the estimators' variance using a Generalized Estimation Equations framework
#'
#' @param mod an object created by \link[corrfuncs]{estimate_model}
#' @param control_arr array of correlation matrices of control group. can be passed as a matrix, each row is a vectorized correlation matrix
#' @param diagnosed_arr same as control_arr but of the diagnosed group
#' @param est_mu whether to use the expected value estimated in the model or the sample average. default TRUE.
#' @return the Fisher information matrix of alpha
#' @seealso \link[corrfuncs]{estimate_model} \link[corrfuncs]{convert_corr_array_to_data_matrix}
#'
#' @export
#'
# todo: compute_sandwich_variance
compute_gee_variance <- function(mod, control_arr, diagnosed_arr, est_mu = TRUE)
  {
  inner <- function(group){
    control_datamatrix <- convert_corr_array_to_data_matrix(control_arr)
    diagnosed_datamatrix <- convert_corr_array_to_data_matrix(diagnosed_arr)
    data <- if(group == 'control') control_datamatrix else diagnosed_datamatrix

    p <- 0.5 + sqrt(1 + 8*ncol(data))/2
    d <- length(mod$alpha)/p

    jacobian <- compute_mu_alpha_jacobian(
      group = group,
      alpha = mod$alpha,
      control_datamatrix = control_datamatrix,
      diagnosed_datamatrix = diagnosed_datamatrix,
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

    solve_sigma <- solve(corrmat_covariance_from_datamatrix(data, est_n = T))
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

  reg_lambda <- mod$regularization$lambda
  reg_p <- mod$regularization$lp

  control_lst <- inner('control')
  diagnosed_lst <- inner('diagnosed')

  I0 <- compute_gee_raw('I0', control_lst) + compute_gee_raw('I0', diagnosed_lst)

  I1 <- compute_gee_raw('I1', control_lst) + compute_gee_raw('I1', diagnosed_lst)
  solve_I0 <- solve(I0)

  out <- solve_I0 %*% I1 %*% solve_I0
  return(out)
}
