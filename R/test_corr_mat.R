test_corr_mat <- function(dta){
  which_not_pos <- which(!apply(with(dta$samples, abind(healthy, sick)), 3, is.positive.definite))
  which_na <- list(
    healthy = which(is.na(dta$healthy), arr.ind = TRUE),
    sick = which(is.na(dta$sick), arr.ind = TRUE)
  )
  return(list(
    which_not_positive_definite = which_not_pos,
    which_na = which_na
  ))
}
