remove_zeros <- function(x, val_to_remove = 0, val_to_replace = 1){
  x[x == val_to_remove] <- val_to_replace
  return(x)
}


matrix_pow <- function(x, pow){
  eig <- eigen(x)
  out <- eig$vectors %*% diag(eig$values ^ pow) %*% t(eig$vectors)
  return(out)
}


mean_sqrt_diag <- function(x) return(mean(sqrt(diag(x))))


sqrt_diag <- function(x) return(sqrt(diag(x)))


force_symmetry <- function(matr) return(0.5 * (t(matr) + matr))


vnorm <- function(x, matr, sqrt = FALSE, solve_matr = FALSE){
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

