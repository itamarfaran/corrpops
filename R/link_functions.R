#' Skeleton of Link Function
#'
#' Skeleton of link function consumed by \link[corrfuncs]{estimate_model}.
#' must have the following elements:
#' - name: optional, a string specifing the link function's name
#' - func: a function with input t (theta), a vectorized correlation matrix, a (alpha) and d, the number of columns in alpha. the output would be a correlation matrix with alpha effect
#' - inverse: a function with input a, that inverts alpha. for example, of we specify func to calculate exp(a), inverse would be log(a)
#' - rev_func: a function with input datamatrix, a vectorized array of correlation matrices, a (alpha) and d, the number of columns in alpha. the output would be the array with the alpha effect reversed
#' - null_value: the null value of the link function
#'
LinkFuncSkeleton <- list(
  name = NA,
  func = function(t, a, d) NA,
  inverse = function(a) NA,
  rev_func = function(datamatrix, a, d) NA,
  null_value = NA
)

#' @describeIn LinkFunctions The multiplicative identity link function
multiplicative_identity <- list(
  name = 'multiplicative_identity',
  func = function(t, a, d) {
    a <- matrix(a, nc = d)
    a_mat <- a %*% t(a)
    diag(a_mat) <- 1
    vector2triangle(t, diag_value = 1)*a_mat
  },
  inverse = function(a) a,
  rev_func = function(datamatrix, a, d){
    a <- matrix(a, nc = d)
    a_mat <- a %*% t(a)
    diag(a_mat) <- 1
    a_vect <- triangle2vector(a_mat)
    return(datamatrix / (rep(1, nrow(datamatrix)) %*% t(a_vect)))
  },
  null_value = 1
)


#' @describeIn LinkFunctions The additive quotent link function
additive_quotent <- list(
  name = 'additive_quotent',
  func = function(t, a, d) {
    a <- as.vector(a)
    a_mat_t <- replicate(length(a), a)
    a_mat <- a_mat_t + t(a_mat_t)
    diag(a_mat) <- 0
    vector2triangle(t, diag_value = 1)/(1 + a_mat)
  },
  inverse = function(a) a,
  rev_func = function(datamatrix, a, d) {
    a <- as.vector(a)
    a_mat_t <- replicate(length(a), a)
    a_mat <- a_mat_t + t(a_mat_t)
    diag(a_mat) <- 0
    a_vect <- triangle2vector(a_mat)
    return(datamatrix * (1 + rep(1, nrow(datamatrix)) %*% t(a_vect)))
  },
  null_value = 0
)


#' Default Link Functions for model
#' @seealso \link[corrfuncs]{LinkFuncSkeleton}
#' @export
LinkFunctions <- list(
  multiplicative_identity = multiplicative_identity,
  additive_quotent = additive_quotent
)
