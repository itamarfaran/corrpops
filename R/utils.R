remove_zeros <- function(x, rem = 0, rep = 1){
  x[x == rem] <- rep
  return(x)
}


matrix_powerer <- function(x, power){
  eig <- eigen(x)
  out <- eig$vectors %*% diag(eig$values ^ power) %*% t(eig$vectors)
  return(out)
}


mean_sqrt_diag <- function(x) return(mean(sqrt(diag(x))))


sqrt_diag <- function(x) return(sqrt(diag(x)))


force_symmetry <- function(matr) return(0.5 * (t(matr) + matr))


vect_norm <- function(x, matr, sqrt = FALSE, solve_matr = FALSE)
{
  if(missing(matr)){
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

