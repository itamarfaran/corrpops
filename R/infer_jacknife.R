infer_jacknife <- function(results)
{
  sick_ind <- as.logical(results$is_sick)

  n_s <- sum(sick_ind)
  n_h <- sum(!sick_ind)

  # estimate_d <- colMeans(results$alpha[sick_ind,])
  const_d <- (n_s - 1)^2 / n_s
  var_d <- const_d * var(results$alpha[sick_ind,])

  # estimate_h <- colMeans(results$alpha[!sick_ind,])
  const_h <- (n_h - 1)^2 / n_h
  var_h <- const_h * var(results$alpha[!sick_ind,])

  # estimate <- (estimate_d + estimate_h)/2
  estimate <- colMeans(results$alpha)
  var_out <- var_d + var_h

  out <- list(
    estimate = estimate,
    variance = var_out
  )
  return(out)
}
