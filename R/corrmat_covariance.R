corrmat_covariance <- function(matr, fisher_z = FALSE, nonpositive = c('stop', 'force', 'ignore'), use_cpp = TRUE){
  nonpositive <- match.arg(nonpositive, c('stop', 'force', 'ignore'))

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
  m <- p * (p-1)/2
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


corrmat_covariance_from_dt <- function(dt, est_n = FALSE, only_diag = TRUE, use_cpp = TRUE, ncores = 1){

  index <- seq_len(nrow(dt))
  func <- function(i) corrmat_covariance(dt[i,], nonpositive = 'ignore', use_cpp = use_cpp)

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
