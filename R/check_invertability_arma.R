check_invertability_arma <- function(coefs, tol = 1e-03)
{
  polyfun <- function(x) 1 - sum(coefs * x ^ (seq_along(coefs)))
  x <- sign(sapply(seq(-1, 1, by = tol), polyfun))
  is_invertable <- all(x == 1) | all(x == -1)
  return(is_invertable)
}
