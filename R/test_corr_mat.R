#' Test if all Matrices of an Array are Positive Definite and without Null
#'
#' @param arr an array of square matrices
#' @return a list of two logicals
#' @export
#'
test_corr_mat <- function(arr){
  out <- list(
    which_not_positive_definite = which(!apply(arr, 3, is.positive.definite)),
    which_na = which(is.na(arr), arr.ind = TRUE)
  )
  return(out)
}
