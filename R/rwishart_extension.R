#' Generate a Covariance Matrix with a Random Effect
#'
#' Generate a Covariance matrix with expected value Sigma, and a random effect controlled by random_effect parameter
#'
#' @param Sigma Expected value of the covariance matrix
#' @param random_effect the control on the random effect size. for 0 < random_effect <= Inf, a covariance matrix with distribution Wishart(Sigma, q)/q is generated with q = p * (1 + 1/random_effect), so a larger random_effect size means a larger variance in the results.
#' @family rWishart2
#' @return a covariance matrix with Expected value Sigma and variance roughly proportional to random_effect
#' @export
generate_random_effect_sigma <- function(Sigma, random_effect = NULL)
{
  if(is.null(random_effect))
    return(Sigma)
  if(random_effect <= 0)
    stop('random_effect must be larger than 0')

  p <- ncol(Sigma)
  p <- p * (1 + 1/random_effect)
  suppressWarnings({out <- rWishart2(1, p, Sigma, random_effect = NULL)[,,1]/p})

  return(out)
}


#' Random Wishart Distributed Matrices with df < p
#'
#' Generate n random matrices, distributed according to the (possibly degenerate) Wishart distribution with parameters Sigma and df, W_p(Sigma, df).
#' with df possibly lower than p, Sigma can be semi positive definite and with a random effect.
#' if df >= p, random_effect is null and Sigma is positive definite, will call \link[stats]{rWishart}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, “degrees of freedom”. can be lower than the dimension of Sigma
#' @param Sigma semi positive definite (p * p) “scale” matrix, the matrix parameter of the distribution.
#' @param random_effect generate a random effect for each matrix with \link[corrfuncs]{generate_random_effect_sigma}
#' @return a numeric array, of dimension p * p * n, where each matrix is semi positive definite covariance matrix, a realization of the (possibly degenerate) Wishart distribution W_p(Sigma, df), where Sigma is possibly an RV itself.
#' @family rWishart2
#' @seealso \link[stats]{rWishart}
#' @export
#'
rWishart2 <- function(n, df, Sigma, random_effect = NULL)
{
  p <- ncol(Sigma)

  if(df >= p & is.null(random_effect) & matrixcalc::is.positive.definite(Sigma))
    return(stats::rWishart(n, df, Sigma))
  if(df < p)
    warning('Wishart degrees of freedom lower than matrix dimension')

  raw_function <- function(k){
    Sigma_ <- if(is.null(Sigma)) Sigma else generate_random_effect_sigma(Sigma, random_effect)
    # todo: not efficient to do an if statement in an inner function
    matrices <- mvtnorm::rmvnorm(n = df, sigma = Sigma_)
    return(t(matrices) %*% matrices)
  }

  out <- lapply(1:n, raw_function)
  out <- simplify2array(out)

  return(out)
}


#' Random Covariance Matrices of Correlated Multivariate Normal RVs
#'
#' Generate n random matrices, distributed proportionally to the (possibly degenerate) Wishart distribution with parameters Sigma and df, W_p(Sigma, df),
#' with df able to be lower than p, Sigma can be semi positive definite and a random effect.
#' Also, unlike the underlying assumption of the Wishart distribution, the covariance matrices are simulated using correlated multivariate Gaussian RVs - not i.i.d.
#' if AR & MA are null, will call \link[corrfuncs]{rWishart2}
#'
#' @inheritParams rWishart2
#' @param AR a vector representing the Auto-regressive part of the multivariate ARMA process
#' @param MA a vector representing the Moving-average part of the multivariate ARMA process
#' @param ncores number of cores to use, if parallelized.
#' @family rWishart2
#' @return a numeric array, of dimension p * p * n, where each matrix is semi positive definite covariance matrix, a realization proportional to the (possibly degenerate) Wishart distribution W_p(Sigma, df), where Sigma is possibly an RV itself.
#' @seealso \link[stats]{rWishart}
#' @export
#'
rWishart_ARMA <- function(n, df, Sigma, random_effect = NULL, AR = NULL, MA = NULL, ncores = 1)
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
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + sum_vector(normal_matrix_arma, (i - arlag):(i - 1), rev(AR[1:arlag]))
        # X[i] = epsilon[i] + X[i-1] + ...
      }
      if(!is.null(MA)){
        malag <- min(max_ma, i - 1)
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + sum_vector(normal_matrix, (i - malag):(i - 1) , rev(MA[1:malag]))
        # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ...
      }
    }
    return(t(normal_matrix_arma) %*% normal_matrix_arma)
  }

  out <- if(ncores > 1) parallel::mclapply(1:n, raw_function, mc.cores = ncores) else lapply(1:n, raw_function)
  out <- simplify2array(out)
  return(out)
}
