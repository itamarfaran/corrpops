#' @describeIn convert_corr_array_to_data_matrix internal function of convert_corr_array_to_data_matrix
convert_corr_array_to_data_matrix_raw <- function(arr){
  matr <- t(apply(arr, 3, triangle2vector))
  return(matr)
}


#' @describeIn convert_data_matrix_to_corr_array internal function of convert_data_matrix_to_corr_array
convert_data_matrix_to_corr_array_raw <- function(matr){
  inner <- function(i) vector2triangle(matr[i,], diag_value = 1)
  arr <- simplify2array(lapply(1:nrow(matr), inner))
  return(arr)
}


#' Convert Array of Correlation Matrices to a Data Matrix
#'
#' Convert an array of Correlation matrices to a data matrix by vectorizing them.
#' Each matrix becomes a row in the output matrix. If arr is already of class matrix, nothing is perfomred.
#'
#' @param arr An array of square, symmetrical matrices or a data matrix of such form.
#' @param verbose should add messages?
#' @return a matrix with the vectorized matrices as rows
#' @family vectriangle
#' @export
convert_corr_array_to_data_matrix <- function(arr, verbose = FALSE)
{
  if(class(arr) == 'array'){
    msg <- 'array transformed from array to matrix'
    out <- convert_corr_array_to_data_matrix_raw(arr)
  } else if (class(arr) %in% c('matrix', 'data.frame')){
    msg <- 'object already in normal data matrix form'
    out <- arr
  } else {
    stop('arr not of class \'array\', \'matrix\' or \'data.frame\'')
  }
  if(verbose) message(msg)
  return(out)
}


#' Convert Data Matrix of Correlation Matrices to an Array
#'
#' Convert a data matrix of vectorized correlation matrices to an array.
#' Each row becomes a matrix in the output array. If matr is already of class array, nothing is performed.
#'
#' @param matr a data matrix of vectorized correlation matrices
#' @param verbose should add messages?
#' @return an array with the rows as correlation matrices
#' @family vectriangle
#' @export
convert_data_matrix_to_corr_array <- function(matr, verbose = FALSE)
{
  if(class(matr) == 'array'){
    msg <- 'object already an array'
    out <- matr
  } else if (class(matr) %in% c('matrix', 'data.frame')){
    msg <- 'matrix transformed from to array'
    out <- convert_data_matrix_to_corr_array_raw(matr)
  } else {
    stop('matr not of class \'array\', \'matrix\' or \'data.frame\'')
  }
  if(verbose) message(msg)
  return(out)
}
