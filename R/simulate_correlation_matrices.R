#' Simulate a Sample of Correlation Matrices
#'
#' will call \link[corrpops]{rWishart_ARMA}. before returning the array, will apply \link[stats]{cov2cor} on each matrix.
#' @inheritParams rWishart_ARMA
#' @return an array of correlation matrices.
#' @export
#'
create_correlation_matrices <- function(n, df, Sigma,
                                        AR = NULL, MA = NULL,
                                        random_effect = NULL,
                                        ncores = 1)
  {
  out <- rWishart_ARMA(n, df, Sigma, AR = AR, MA = MA, random_effect = random_effect, ncores = ncores)
  out <- lapply(1:n, function(b) force_symmetry(stats::cov2cor(out[,,b])))
  out <- simplify2array(out)
  return(out)
}
