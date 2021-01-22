convert_corr_array_to_data_matrix_raw <- function(arr) t(apply(arr, 3, triangle2vector))


convert_corr_array_to_data_matrix <- function(obj, verbose = FALSE)
{
  if(class(obj) == 'array'){
    msg <- 'obj transformed from array to matrix'
    out <- convert_corr_array_to_data_matrix_raw(obj)
  } else if (class(obj) %in% c('matrix', 'data.frame')){
    msg <- 'obj already in normal data matrix form'
    out <- obj
  } else {
    stop('obj not of class \'array\', \'matrix\' or \'data.frame\'')
  }
  if(verbose) message(msg)
  return(out)
}
