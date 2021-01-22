test_corr_mat <- function(dta){
  out <- list(
    which_not_positive_definite = which(!apply(dta, 3, is.positive.definite)),
    which_na = which(is.na(dta), arr.ind = TRUE)
  )
  return(out)
}
