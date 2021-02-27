#' Check Invertibility/Stationarity of an AR/MA Process
#'
#'  Check if an AR process is stationary or an MA process is invertible,
#'  by testing if the polynomial defined by it's coefficients doesn't have a root between -1 and 1
#'
#' @param coefs a numeric vector with the AR (or MA) coefficients
#' @param tol the tolerance search for the roots.
#' @return true if the process is invertible/stationary.
#'
#' @export
is_invertable_arma <- function(coefs, tol = 1e-03)
{
  polyfun <- function(x) 1 - sum(coefs * x ^ (seq_along(coefs)))
  x <- sign(sapply(seq(-1, 1, by = tol), polyfun))
  is_invertable <- all(x == 1) | all(x == -1)
  return(is_invertable)
}
