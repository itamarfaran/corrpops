test_that('sum_vector of ones gets the sum', {
  n <- 100
  p <- 10
  matr <- matrix(1, n, p)
  expect_equal(sum_vector(matr), rep(n, p))
  expect_equal(sum_vector(matr, by_row = FALSE), rep(p, n))
})

test_that('calculate_mean_matrix calculates the mean matrix', {
  arr <- array(dim = c(3, 3, 3))
  arr[,,1] <- 1
  arr[,,2] <- 2
  arr[,,3] <- 3

  expect_equal(calculate_mean_matrix(arr), matrix(2, 3, 3))
})
