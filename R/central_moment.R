central_moment <- function(x, norm = TRUE) {
  out <- numeric(4)
  names(out) <- c("mean", "variance", "skewness", ifelse(norm, "kurtosis", "exkurtosis"))

  n <- length(x)
  out[1] <- m <- mean(x)
  out[2] <- v <- var(x)

  s <- sqrt(v * (1 - 1/n))
  norm_x <- (x - m)/s

  skew <- mean(norm_x ^ 3)
  kurt <- mean(norm_x ^ 4)

  out[3] <- (sqrt(n*(n - 1)) / (n - 2)) * skew
  out[4] <- ( (n - 1) / ((n - 2) * (n - 3)) ) * ((n + 1) * kurt + 6) - 3*norm

  return(out)
}
