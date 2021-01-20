vector_sum <- function(matr, index, constants, weights = FALSE, by_row = TRUE)
  {
  if(!by_row) matr <- t(matr)

  out <- numeric(ncol(matr))
  if(missing(index)) index <- 1:nrow(matr)
  if(missing(constants)) constants <- rep(1, length(index))
  if(weights) constants <- constants / sum(constants)

  for(i in seq_along(index)) out <- out + constants[i]*matr[index[i],]
  return(out)
}


matrix_sum <- function(array_, index, constants, weights = FALSE)
  {
  dim_array <- dim(array_)
  out <- matrix(0, nrow = dim_array[1], ncol = dim_array[2])

  if(missing(index)) index <- 1:dim_array[3]
  if(missing(constants)) index <- rep(1, length(index))
  if(weights) constants <- constants/sum(constants)

  for(i in seq_along(index)) out <- out + constants[i]*array_[,, index[i]]

  return(out)
}


calculate_mean_matrix <- function(matrix_array) matrix_sum(matrix_array, weights = TRUE)

