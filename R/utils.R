#' Replace a Value in a vector/array
#'
#' @param x the array
#' @param rem the value to replace, default 0
#' @param rep the value to replace with, default 1
#' @return  the array with rem replaces with rep
#' @export
remove_zeros <- function(x, rem = 0, rep = 1){
  x[x == rem] <- rep
  return(x)
}


#' Raise a Square Matrix to a Non-Integer Power
#'
#' Raise a Square Matrix to a non-integer power using it's eigen decomposition
#' @param x a square matrix
#' @param power a real number
#' @return the matrix raised to the power of power
#' @export
matrix_power <- function(x, power){
  eig <- eigen(x)
  out <- eig$vectors %*% diag(eig$values ^ power) %*% t(eig$vectors)
  return(out)
}


#' Get the Average of the Square-rooted Diagonal of a Matrix
#' @param x a square matrix
#' @export
mean_sqrt_diag <- function(x) return(mean(sqrt(diag(x))))


#' Get the Square-rooted Diagonal of a Matrix
#' @param x a square matrix
#' @export
sqrt_diag <- function(x) return(sqrt(diag(x)))


#' Force a Matrix to be Symmetrical
#'
#' Force a matrix to be symmetrical by averaging between the matrix and it's transpose
#' @export
force_symmetry <- function(matr) return(0.5 * (t(matr) + matr))


#' Calculate the Mahalonobis Norm of a Vector
#' @param x a vector
#' @param matr the weighting matrix in the Mahalonobis norm. if null, the l2 norm is calculated.
#' @param sqrt if false, don't take the sqrt of the product
#' @param solve_matr if true, invert matr before calculating the norm
#' @return the norm, a real number >= 0
#' @export
vect_norm <- function(x, matr = NULL, sqrt = TRUE, solve_matr = FALSE)
{
  if(is.null(matr)){
    out <- sum(x^2)
  } else {
    if(solve_matr)
      matr <- solve(matr)
    out <- as.vector(t(x) %*% matr %*% x)
  }

  if(sqrt)
    out <- sqrt(out)

  return(out)
}

