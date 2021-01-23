multiplicative_identity <- list(
  name = 'multiplicative_identity',
  func = function(t, a, d) {
    a <- matrix(a, nc = d)
    a_mat <- a %*% t(a)
    diag(a_mat) <- 1
    vector2triangle(t, diag_value = 1)*a_mat
  },
  inverse = function(a) a,
  rev_func = function(dt, a, d){
    a <- matrix(a, nc = d)
    a_mat <- a %*% t(a)
    diag(a_mat) <- 1
    a_vect <- triangle2vector(a_mat)
    return(dt / (rep(1, nrow(dt)) %*% t(a_vect)))
  },
  null_value = 1
)


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
  rev_func = function(dt, a, d) {
    a <- as.vector(a)
    a_mat_t <- replicate(length(a), a)
    a_mat <- a_mat_t + t(a_mat_t)
    diag(a_mat) <- 0
    a_vect <- triangle2vector(a_mat)
    return(dt * (1 + rep(1, nrow(dt)) %*% t(a_vect)))
  },
  null_value = 0
)


LinkFunctions <- list(
  multiplicative_identity = multiplicative_identity,
  additive_quotent = additive_quotent
)
