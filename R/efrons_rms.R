efron_bias_correction <- function(rms, p){
  rms <- sqrt(p/(p - 1) * (rms ^ 2 - 1/(p - 1)))
  return(rms)
}


efrons_rms_sample <- function(df, p = NULL){
  df <- convert_corr_array_to_data_matrix(df)
  rms_mean <- mean(sqrt(rowMeans(df^2)))

  if(!is.null(p))
    rms_mean <- efron_bias_correction(rms_mean)

  return(rms_mean)
}


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


efrons_effective_sample_size <- function(n, rms){
  out <- n/(1 + (n - 1) * rms^2)
  return(out)
}
