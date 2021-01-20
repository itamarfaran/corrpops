regularize_matrix <- function(matr,
                              const = 1,
                              method = c("constant", "avg_diag", "increase_diag"),
                              only_if_singular = TRUE)
  {

  method <- match.arg(method, c("constant", "avg_diag", "increase_diag"))
  if(nrow(matr) != ncol(matr))
    stop("Matrix is not square.")

  if(only_if_singular)
    if(!matrixcalc::is.singular.matrix(matr))
      return(matr)

  p <- nrow(matr)
  if(method == "constant"){
    if(const < 0) stop("in method 'constant' const must be greater or equal to 0")
    out <- matr + const * diag(p)
  }
  if(method == "avg_diag"){
    if(const < 0 | const > 1) stop("in method 'avg_diag' const must be between 0-1")
    out <- (1 - const) * matr + const * mean(diag(matr)) * diag(p)
  }
  if(method == "increase_diag"){
    if(const < 0 | const > 1) stop("in method 'increase_diag' const must be between 0-1")
    out <- (1 - const) * matr + const * diag(diag(matr))
  }
  return(out)
}
