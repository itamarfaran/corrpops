#' Regularize Matrix
#'
#' Regularize a square matrix by increasing it's diagonal
#' @param matr the matrix to be regularized
#' @param const the constant
#' @param only_if_singular if true, will regularize only if the matrix is singular
#' @return if method = 'constant', will add the constant to the diagonal. if method = 'avg_diag', will calculate a weighted average of the original matrix and an identity matrix multiplied by the avregae value of the diagonal, with weights 1 - const and const. if method = 'avg_diag', will calculate a weighted average of the original matrix and the diagonal of the matrix, with weights 1 - const and const
#'
#' @export
#'
regularize_matrix <- function(matr, const = 1,
                              method = c('constant', 'avg_diag', 'increase_diag'),
                              only_if_singular = TRUE)
{

  method <- match.arg(method, c('constant', 'avg_diag', 'increase_diag'))
  if(nrow(matr) != ncol(matr))
    stop('matrix is not square')

  if(only_if_singular)
    if(!matrixcalc::is.singular.matrix(matr))
      return(matr)

  p <- nrow(matr)
  if(method == 'constant'){
    if(const < 0) stop('in method \'constant\' const must be greater or equal to 0')
    out <- matr + const * diag(p)
  }
  if(method == 'avg_diag'){
    if(const < 0 | const > 1) stop('in method \'avg_diag\' const must be between 0-1')
    out <- (1 - const) * matr + const * mean(diag(matr)) * diag(p)
  }
  if(method == 'increase_diag'){
    if(const < 0 | const > 1) stop('in method \'increase_diag\' const must be between 0-1')
    out <- (1 - const) * matr + const * diag(diag(matr))
  }
  return(out)
}
