#' Bias Correction of Efron's RMS
#' @param rms the RMS factor
#' @param p the dimension of the matrix
#' @family efron_rms
#' @export
efron_bias_correction <- function(rms, p){
  rms <- sqrt(p/(p - 1) * (rms ^ 2 - 1/(p - 1)))
  return(rms)
}


#' Calculate the Average Efron's RMS
#' Calculate the average Efron's RMS over a sample of corraltion matrices
#' @param arr an array of correlation matrices (can be in data matrix-form)
#' @param p the dimension of the matrix
#' @family efron_rms
#' @export
efrons_rms_sample <- function(arr, p = NULL){
  arr <- convert_corr_array_to_data_matrix(arr)
  rms_mean <- mean(sqrt(rowMeans(arr^2)))

  if(!is.null(p))
    rms_mean <- efron_bias_correction(rms_mean)

  return(rms_mean)
}


#' Calculate Efron's RMS Over A Correlation Matrix
#' @param m a correlation matrix
#' @param p the dimension of the matrix. defult null. if not null, conduct a bias correction.
#' @references Bradley Efron. Large-scale inference: empirical Bayes methods for estimation, testing, and prediction, volume 1. Cambridge University Press, 2012.
#' @family efron_rms
#' @export
efrons_rms <- function(m, p = NULL){
  if(ncol(m) != nrow(m))
    stop('m is not square')
  if(any(diag(m) != 1))
    stop('diag of m is not 1')

  m_vect <- triangle2vector(m, diag = FALSE)
  rms <- sqrt(mean(m_vect^2))

  if(!is.null(p))
    rms <- efron_bias_correction(rms)
  return(rms)
}


#' Calculate Efron's Effective Sample Size
#'
#' Estimate the Effective Sample Size, reduced due to cross-row correlation.
#' @param n the number of ovservations
#' @param rms Efron's Root Mean Square coefficient.
#' @family efron_rms
#' @export
efrons_effective_sample_size <- function(n, rms){
  out <- n/(1 + (n - 1) * rms^2)
  return(out)
}
