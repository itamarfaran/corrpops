#' Check Invertibility/Stationarity of an AR/MA Process
#'
#'  Check if an AR proccess is stationary, or an MA process is invertable,
#'  by ...
#'
#' @param coefs the AR (or MA) coefficients. numeric vector
#' @param tol the tolerance to check the ...
#' @return logical
#'
#' @export
is_invertable_arma <- function(coefs, tol = 1e-03)
{
  polyfun <- function(x) 1 - sum(coefs * x ^ (seq_along(coefs)))
  x <- sign(sapply(seq(-1, 1, by = tol), polyfun))
  is_invertable <- all(x == 1) | all(x == -1)
  return(is_invertable)
}
