infer_jacknife <- function(results){
  sick_ind <- as.logical(results$is_sick)

  n_s <- sum(sick_ind)
  n_h <- sum(!sick_ind)

  estimate_d <- colMeans(results$alpha[sick_ind,])
  const_d <- (n_s - 1)^2/n_s
  var_d <- var(results$alpha[sick_ind,])*const_d

  estimate_h <- colMeans(results$alpha[!sick_ind,])
  const_h <- (n_h - 1)^2/n_h
  var_h <- var(results$alpha[!sick_ind,])*const_h

  estimate <- colMeans(results$alpha)
  # estimate <- (estimate_d + estimate_h)/2
  var_out <- var_d + var_h

  return(list(
    estimate = estimate,
    variance = var_out
  ))
}
