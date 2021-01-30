#' @describeIn convert_corr_array_to_data_matrix internal function of convert_corr_array_to_data_matrix
convert_corr_array_to_data_matrix_raw <- function(arr) t(apply(arr, 3, triangle2vector))


#' Convert Array of Symmetrical Matrices to a Data Matrix
#'
#' Convert an array of symmetrical matrices to a data matrix by vectorizing them.
#' Each matrix becomes a row in the output matrix. If obj is already of class matrix, nothing is perfomred.
#' @seealso \link[corrfuncs]{triangle2vector}
#'
#' @param obj An array of square, symmetrical matrices or a data matrix of such form.
#' @return Matrix with vectorized matrices as rows
#' @export
convert_corr_array_to_data_matrix <- function(obj, verbose = FALSE)
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
