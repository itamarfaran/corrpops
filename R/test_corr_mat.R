#' Test if all Matrices of an Array are Positive Definite and without Missing Values
#'
#' @param arr an array of square matrices
#' @return a list of index vectors:
#' @return - which_not_positive_definite: which matrices in the array aren't positive_definite
#' @return - which_na: which matrices in the array have missing values
#' @export
#'
test_corr_mat <- function(arr){
  out <- list(
    which_not_positive_definite = which(!apply(arr, 3, matrixcalc::is.positive.definite)),
    which_na = which(is.na(arr), arr.ind = TRUE)
  )
  return(out)
}
