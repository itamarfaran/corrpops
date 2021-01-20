create_correlation_matrices <- function(real_corr,
                                        real_var,
                                        sample_size, df = 0,
                                        AR = NULL, MA = NULL,
                                        random_effect = NULL,
                                        ncores = 1)
  {
  Sigma <- real_corr
  out <- rWishart_ARMA(sample_size, df, Sigma, AR = AR, MA = MA, random_effect = random_effect, ncores = ncores)
  out <- lapply(1:sample_size, function(b) force_symmetry(cov2cor(out[,,b])))
  out <- simplify2array(out)
  return(out)
}
