#' Transform a Symmetric Matrix to it's Vectorised Form and Vise Versa
#'
#' A symmetric matrix can be entirely summarized with one of its triangles.
#' these functions transform a matrix's lower triangle to a vectorized form, and transform a vector to a symmetric matrix.
#' In correlation matrices, the diagonal is redundant, so it is possible to vectorize the off diagonal triangle.
#'
#' @param matr the matrix to vectorize it's lower triangle
#' @param vect the vector to transform to a symmetric matrix
#' @param diag if true, vectorize (or construct) the diagonal of the matrix as well
#' @param diag_value if diag is set to false, what value to place on the diagonal. default is NA, for correlation matrices set to 1.
#' @return the transformed matrix / vector
#' @family vectriangle
#' @export
#'
triangle2vector <- function(matr, diag = FALSE){
  if(nrow(matr) != ncol(matr))
    stop('matrix is not square')

  out <- as.vector(matr[lower.tri(matr, diag = diag)])

  return(out)
}


#' @describeIn triangle2vector inverse of triangle2vector
#' @family vectriangle
#' @export
vector2triangle <- function(vect, diag = FALSE, diag_value = NA){
  m <- length(vect)

  one <- ifelse(diag, -1, 1)
  p <- 0.5 * c(one + sqrt(1 + 8 * m), one - sqrt(1 + 8 * m))
  p <- p[which((p == round(p)) & p == abs(p))]

  if(length(p) == 0)
    stop('vect length does not fit size of triangular matrix')

  out <- matrix(0, ncol = p, nrow = p)
  out[lower.tri(out, diag = diag)] <- vect

  out <- out + t(out) - diag(diag(out))
  if(!diag) diag(out) <- diag_value

  return(out)
}
