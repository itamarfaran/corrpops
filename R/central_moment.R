#' Calculate Central Moments
#'
#' Calculate the sample Mean, Variance, Skewness and Kurtosis
#'
#' @param x numeric vector
#' @param exkurtosis logical, whther to caculate the ex-kurtosis or kurtosis. default True
#' @return named numeric vector of length 4
#'
#' @export
#'
central_moments <- function(x, exkurtosis = TRUE)
{
  out <- numeric(4)
  names(out) <- c('mean', 'variance', 'skewness', ifelse(exkurtosis, 'kurtosis', 'exkurtosis'))

  n <- length(x)
  out[1] <- m <- mean(x)
  out[2] <- v <- var(x)
  s <- sqrt(v * (1 - 1/n))
  x_normalized <- (x - m)/s

  skew <- mean(x_normalized ^ 3)
  kurt <- mean(x_normalized ^ 4)

  out[3] <- (sqrt(n*(n - 1)) / (n - 2)) * skew
  out[4] <- ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * kurt + 6) - 3 * exkurtosis

  return(out)
}
