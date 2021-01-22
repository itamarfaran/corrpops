triangle2vector <- function(matr, diag = FALSE){
  if(nrow(matr) != ncol(matr))
    stop('matrix is not square')

  out <- as.vector(matr[lower.tri(matr, diag = diag)])

  return(out)
}


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
