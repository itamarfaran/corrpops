#' Simulate a Sample of Correlation Matrices
#'
#'
#' @inheritParams rWishart_ARMA
#' @return an array of correlation matrices.
#' @export
#'
create_correlation_matrices <- function(n, df, Sigma,
                                        AR = NULL, MA = NULL,
                                        random_effect = NULL,
                                        ncores = 1)
  {
  Sigma <- real_corr
  out <- rWishart_ARMA(n, df, Sigma, AR = AR, MA = MA, random_effect = random_effect, ncores = ncores)
  out <- lapply(1:sample_size, function(b) force_symmetry(cov2cor(out[,,b])))
  out <- simplify2array(out)
  return(out)
}
