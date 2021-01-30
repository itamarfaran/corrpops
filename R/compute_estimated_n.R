#' @describeIn compute_estimated_n internal function of compute_estimated_n
compute_estimated_n_raw <- function(est, theo, only_diag = FALSE)
{
  if(only_diag){
    x <- diag(theo)
    y <- diag(est)
  } else {
    x <- triangle2vector(theo, diag = TRUE)
    y <- triangle2vector(est, diag = TRUE)
  }
  out <- lm(x ~ 0 + y)$coef
  return(out)
}

#' Estimate Degrees of Freedom by Linear Projection
#'
#' Estimate the Degrees of Freedom of a model by projecting the theoretical variance
#' on the estimated (empirical) variance.
#' When using compute_estimated_n, the theoretical variance is computed using corrmat_covariance
#'
#' @param obj array of correlation matrices or a matrix of diagnolized correlation matrices.
#' @param est the empirical covariance matrix, \eqn{n^{-1}\left(Y-\hat{\mu}\right)\left(Y-\hat{\mu}\right)^{t}}
#' @param theo the theoretical covariance matrix
#' @param only_diag logical, whether to use only the diagonal (variances) or the whole matrix. Default FALSE
#' @return numeric, the estimated degrees of freedom
compute_estimated_n <- function(obj, only_diag = TRUE)
{
  datamatrix <- convert_corr_array_to_data_matrix(obj)
  est <- t(datamatrix) %*% datamatrix / nrow(datamatrix)
  theo <- corrmat_covariance_from_datamatrix(datamatrix)
  out <- compute_estimated_n_raw(est = est, theo = theo, only_diag = only_diag)
  return(out)
}
