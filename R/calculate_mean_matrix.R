#' Sum Over Matrix Rows / Columns
#'
#' Sum over matrix rows / columns
#'
#' @param matr matrix to sum over rows / columns
#' @param index index of rows / columns to sum. default is to sum over all rows.
#' @param constants constants to multiply the indices before summation. default is 1 for all constants
#' @param weights should treat constants as weights (sum of all constants = 1)? default false
#' @param by_row should sum over rows (default) or by columns?
#' @return a vector
#'
#' @export
#'
sum_vector <- function(matr, index, constants, weights = FALSE, by_row = TRUE)
  {
  if(!by_row)
    matr <- t(matr)

  out <- numeric(ncol(matr))

  if(missing(index))
    index <- 1:nrow(matr)
  if(missing(constants))
    constants <- rep(1, length(index))
  if(weights)
    constants <- constants / sum(constants)

  for(i in seq_along(index))
    out <- out + constants[i] * matr[index[i],]

  return(out)
}


#' Sum Over Array's Matrices
#'
#' Sum pver array's matrices
#'
#' @param arr array to sum over 3rd axis.
#' @param index index of matrices to sum. default is to sum over all matrices.
#' @param constants constants to multiply the matrices before summation. default is 1 for all constants
#' @param weights should treat constants as weights (sum of all constants = 1)? default false
#' @return a matrix
#'
#' @export
#'
sum_matrix <- function(arr, index, constants, weights = FALSE)
  {
  dim_arr <- dim(arr)
  out <- matrix(0, nrow = dim_arr[1], ncol = dim_arr[2])

  if(missing(index))
    index <- 1:dim_arr[3]
  if(missing(constants))
    index <- rep(1, length(index))
  if(weights)
    constants <- constants/sum(constants)

  for(i in seq_along(index))
    out <- out + constants[i] * arr[,, index[i]]

  return(out)
}


#' Calculate Mean Matrix of an Array
#'
#' Calculate mean matrix of an array
#'
#' @param arr array to sum over 3rd axis.
#' @return averaged matrix
#'
#' @export
#'
calculate_mean_matrix <- function(arr){
  return(sum_matrix(arr, weights = TRUE))
}
