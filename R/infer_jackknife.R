#' Infer Jackknife Results
#'
#' Calculate the Jackknife Variance based on Jackknife Estimation
#'
#' @param obj output object of \link[corrpops]{estimate_model_jacknife}
#' @return a list consisting of estimates and variance matrix
#' @return - estimate: Jacknife estimate of alpha (average estimate over all jackknife iterations)
#' @return - variance: Jacknife estimate of the variance matrix
#'
#' @export

infer_jackknife <- function(obj)
{
  diagnosed_ind <- as.logical(obj$is_diagnosed)

  n_s <- sum(diagnosed_ind)
  n_h <- sum(!diagnosed_ind)

  # estimate_d <- colMeans(obj$alpha[diagnosed_ind,])
  const_d <- (n_s - 1)^2 / n_s
  var_d <- const_d * stats::var(obj$alpha[diagnosed_ind,])

  # estimate_h <- colMeans(obj$alpha[!diagnosed_ind,])
  const_h <- (n_h - 1)^2 / n_h
  var_h <- const_h * stats::var(obj$alpha[!diagnosed_ind,])

  # estimate <- (estimate_d + estimate_h)/2
  estimate <- colMeans(obj$alpha)
  var_out <- var_d + var_h

  out <- list(
    estimate = estimate,
    variance = var_out
  )
  return(out)
}
