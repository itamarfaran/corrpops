#' Calculate the Covariance Matrix of a Correlation Matrix
#'
#' Calculate the covariance matrix of a correlation matrix. The matrix can be in vectorized form.
#'
#' @param matr the correlation matrix. can be vectorized.
#' @param fisher_z if true, calculate the covariance matrix of a fisher-Z transformed correlation matrix. It is assumed that the correlations are Fisher transformed, and the matrix will under go the Inverse Fisher Transformation.
#' @param nonpositive error handling when matrix with not positive definite. can be one of 'stop', 'force', 'ignore'. if 'force' is chosen, \link[Matrix]{nearPD} will be used.
#' @param use_cpp, whether to use c++ source code. default true.
#' @return the covariance matrix
#' @seealso \link[corrfuncs]{corrmat_covariance_from_dt}
#' @export

corrmat_covariance <- function(matr, fisher_z = FALSE, nonpositive = c('stop', 'force', 'ignore'), use_cpp = TRUE){
  nonpositive <- match.arg(nonpositive, c('stop', 'force', 'ignore'))

  if(fisher_z)
    matr <- (exp(2*matr) - 1)/(exp(2*matr) + 1)

  if(is.vector(matr))
    matr <- vector2triangle(matr, diag_value = 1)

  if(nonpositive != 'ignore'){
    if(!matrixcalc::is.positive.definite(matr)){
      if(nonpositive == 'force') {
        matr <- as.matrix(Matrix::nearPD(matr, corr = TRUE, doSym = TRUE)$mat)
        warning('forced positive definiteness with Matrix::nearPD')
      } else{
        stop('matr not positive definite')
      }
    }
  }

  p <- nrow(matr)
  m <- 0.5 * p * (p - 1)
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i))) - use_cpp
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p)) - use_cpp

  corrcalc <- if(use_cpp) corrcalc_c else corrcalc_r
  output <- corrcalc(matr, m, order_vecti, order_vectj)
  output <- output + t(output) - diag(diag(output))

  if(fisher_z){
    diagonalized <- triangle2vector(matr)
    transformed <- 1/(1 - diagonalized ^ 2)
    gradient <- diag(transformed)
    output <- gradient %*% output %*% gradient
  }
  return(output)
}


#' Calculate the Average Covariance Matrix of a Sample of Correlation Matrix
#'
#' Calculate the average covariance matrix of multiple correlation matrices.
#'
#' @param dt the sample of correlation matrices. can be an array of matrices, or a matrix of vectorized correlation matrices.
#' @param fisher_z if true, calculate the covariance matrix of a fisher-Z transformed correlation matrix. It is assumed that the correlations are Fisher transformed, and the matrix will under go the Inverse Fisher Transformation.
#' @param est_n whether to divide the covariance matrix by the estimated degrees of freedom, with a linear projection.
#' @param only_diag passed to \link[corrfuncs]{compute_estimated_n}
#' @param nonpositive error handling when matrix with not positive definite. can be one of 'stop', 'force', 'ignore'. if 'force' is chosen, \link[Matrix]{nearPD} will be used.
#' @param use_cpp whether to use c++ source code. default true.
#' @param ncores number of cores to use, if parallelization is in use.
#' @return the averaged covariance matrix
#' @seealso \link[corrfuncs]{corrmat_covariance}
#' @export
corrmat_covariance_from_dt <- function(dt, fisher_z = FALSE,
                                       est_n = FALSE, only_diag = TRUE,
                                       nonpositive = c('ignore', 'stop', 'force'),
                                       use_cpp = TRUE, ncores = 1){

  dt <- convert_corr_array_to_data_matrix(dt)
  index <- seq_len(nrow(dt))
  func <- function(i) corrmat_covariance(dt[i,], fisher_z = fisher_z, nonpositive = nonpositive, use_cpp = use_cpp)

  if(ncores > 1){
    covariances <- parallel::mclapply(index, func, mc.cores = ncores)
  } else {
    covariances <- lapply(index, func)
  }

  sigma <- calculate_mean_matrix(simplify2array(covariances))

  if(est_n){
    est <- t(dt) %*% dt / nrow(dt)
    estimated_n <- compute_estimated_n_raw(est = est, theo = sigma, only_diag = only_diag)
    sigma <- sigma/estimated_n
  }
  return(sigma)
}
