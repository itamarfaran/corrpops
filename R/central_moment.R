#' Calculate Central Moments
#'
#' Calculate the unbiased estimates of a (univariate) sample's Mean, Variance, Skewness and Kurtosis
#'
#' @param x a numeric vector
#' @param exkurtosis if true, calculates the ex-kurtosis instead of the kurtosis (subtract 3 from the kurtosis, to ensure the kurtosis of a Gaussian RV is 0). default True
#' @return a named numeric vector consisting of the estimates for the first four centralized moments.
#'
#' @export
#'
central_moments <- function(x, exkurtosis = TRUE)
{
  out <- numeric(4)
  names(out) <- c('mean', 'variance', 'skewness', ifelse(exkurtosis, 'kurtosis', 'exkurtosis'))

  n <- length(x)
  out[1] <- m <- mean(x)
  out[2] <- v <- stats::var(x)
  s <- sqrt(v * (1 - 1/n))
  x_normalized <- (x - m)/s

  skew <- mean(x_normalized ^ 3)
  kurt <- mean(x_normalized ^ 4)

  out[3] <- (sqrt(n*(n - 1)) / (n - 2)) * skew
  out[4] <- ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * kurt + 6) - 3 * exkurtosis

  return(out)
}
