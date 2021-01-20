convert_corr_array_to_data_matrix <- function(array_) t(apply(array_, 3, triangle2vector))


convert_corr_array_to_data_matrix_test <- function(obj, verbose = FALSE){
  if(class(obj) == "array"){
    message_ <- 'obj transformed from array to matrix'
    out <- convert_corr_array_to_data_matrix(obj)
  } else if (class(obj) %in% c("matrix", "data.frame")){
    message_ <- 'obj already in normal data matrix form'
    out <- obj
  } else {
    stop('obj not of class "array", "matrix" or "data.frame"')
  }
  if(verbose) message(message_)
  return(out)
}
