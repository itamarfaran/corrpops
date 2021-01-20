check_invertability_arma <- function(coefs, perc = 1e-03){
  polyfun <- function(x) 1 - sum(coefs * x ^ (seq_along(coefs)))
  x <- sign(sapply(seq(-1, 1, by = perc), polyfun))
  return(all(x == 1) | all(x == -1))
}
