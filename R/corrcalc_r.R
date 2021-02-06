#' Calculate the Covariance Matrix of a Vectorized Correlation matrix
#'
#' For inner use.
#'
#' @seealso \link[corrfuncs]{corrmat_covariance}
#' @seealso \link[corrfuncs]{triangle2vector}
corrcalc_r <- function(matr, m, order_vecti, order_vectj)
{
  output <- matrix(0, nrow = m, ncol = m)

  for(i1 in 1:m){
    for(j1 in i1:m){
      i <- order_vecti[i1]
      j <- order_vectj[i1]
      k <- order_vecti[j1]
      l <- order_vectj[j1]

      matr_ij <- matr[i,j]
      matr_kl <- matr[k,l]
      matr_ik <- matr[i,k]
      matr_il <- matr[i,l]
      matr_jk <- matr[j,k]
      matr_jl <- matr[j,l]

      output[i1,j1] <-
        (matr_ij*matr_kl/2) * (matr_ik^2 + matr_il^2 + matr_jk^2 + matr_jl^2) -
        matr_ij*(matr_ik*matr_il + matr_jk*matr_jl) -
        matr_kl*(matr_ik*matr_jk + matr_il*matr_jl) +
        (matr_ik*matr_jl + matr_il*matr_jk)
    }
  }
  return(output)
}
