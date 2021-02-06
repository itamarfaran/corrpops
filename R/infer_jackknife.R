#' Infer Jackknife Results
#'
#' Claculate the Jacknife Variance base on Jackknife Estimation
#'
#' @param results output object of \link[corrfuncs]{estimate_model_jacknife}
#' @return list of estimates and variance matrix
#'
#' @export

infer_jackknife <- function(results)
{
  diagnosed_ind <- as.logical(results$is_diagnosed)

  n_s <- sum(diagnosed_ind)
  n_h <- sum(!diagnosed_ind)

  # estimate_d <- colMeans(results$alpha[diagnosed_ind,])
  const_d <- (n_s - 1)^2 / n_s
  var_d <- const_d * var(results$alpha[diagnosed_ind,])

  # estimate_h <- colMeans(results$alpha[!diagnosed_ind,])
  const_h <- (n_h - 1)^2 / n_h
  var_h <- const_h * var(results$alpha[!diagnosed_ind,])

  # estimate <- (estimate_d + estimate_h)/2
  estimate <- colMeans(results$alpha)
  var_out <- var_d + var_h

  out <- list(
    estimate = estimate,
    variance = var_out
  )
  return(out)
}
