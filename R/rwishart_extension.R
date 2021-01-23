generate_random_effect_sigma <- function(Sigma, random_effect = NULL)
{
  if(is.null(random_effect))
    return(Sigma)

  p <- ncol(Sigma)
  if(0 <= random_effect & random_effect < Inf)
    p <- p * (1 + 1/random_effect)

  out <- rWishart(1, p, Sigma)[,,1]/p

  return(out)
}


rWishart2 <- function(n = 1, df, Sigma, random_effect = NULL)
{
  p <- ncol(Sigma)

  if(df >= p & is.null(random_effect))
    return(rWishart(n, df, Sigma))
  if(df < p)
    warning('Wishart degrees of freedom lower than matrix dimension')

  raw_function <- function(k)
  {
    Sigma <- generate_random_effect_sigma(Sigma, random_effect)
    matrices <- mvtnorm::rmvnorm(n = df, sigma = Sigma)
    return(t(matrices) %*% matrices)
  }

  out <- lapply(1:n, raw_function)
  out <- simplify2array(out)

  return(out)
}


rWishart_ARMA <- function(n = 1, df, Sigma, AR = NULL, MA = NULL, random_effect = NULL, ncores = 1)
{
  p <- ncol(Sigma)

  if(is.null(MA) & is.null(AR))
    return(rWishart2(n = n, df = df, Sigma = Sigma, random_effect = random_effect))
  if(df < p)
    warning('Wishart degrees of freedom lower than matrix dimension')

  max_ar <- max_ma <- NULL

  if(!is.null(AR)){
    if(!is_invertable_arma(AR))
      stop('AR process not stationary')
    max_ar <- length(AR)
  }
  if(!is.null(MA)){
    if(!is_invertable_arma(MA))
      stop('MA process not invertable')
    max_ma <- length(MA)
  }

  raw_function <- function(k){
    Sigma <- generate_random_effect_sigma(Sigma, random_effect)
    normal_matrix_arma <- mvtnorm::rmvnorm(n = df, sigma = Sigma)
    normal_matrix <- mvtnorm::rmvnorm(n = df, sigma = Sigma)

    for(i in 2:df){
      if(!is.null(AR)){
        arlag <- min(max_ar, i - 1)
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + vector_sum(normal_matrix_arma, (i - arlag):(i - 1), rev(AR[1:arlag]))
        # X[i] = epsilon[i] + X[i-1] + ...
      }
      if(!is.null(MA)){
        malag <- min(max_ma, i - 1)
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + vector_sum(normal_matrix, (i - malag):(i - 1) , rev(MA[1:malag]))
        # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ...
      }
    }
    return(t(normal_matrix_arma) %*% normal_matrix_arma)
  }

  out <- if(ncores > 1) parallel::mclapply(1:n, raw_function, mc.cores = ncores) else lapply(1:n, raw_function)
  out <- simplify2array(out)
  return(out)
}
