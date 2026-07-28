# Dense-block and structured Schur sandwich estimators for corrpops
#
# These functions are for the one-dimensional multiplicative-identity model
#
#   mu_H(theta, alpha) = theta,
#   mu_D(theta, alpha) = theta * e(alpha),
#   e_jk(alpha)        = alpha_j alpha_k.
#
# The two public functions use the same joint estimating equations. The dense
# version forms the complete (m + p) sandwich; the Schur version eliminates the
# m-dimensional theta block before forming the sandwich.


# -----------------------------------------------------------------------------
# Shared algebra
# -----------------------------------------------------------------------------

.block_effect_terms <- function(alpha, pair_index) {
  m <- nrow(pair_index)
  p <- length(alpha)
  edge <- seq_len(m)
  j <- pair_index[, 1L]
  k <- pair_index[, 2L]

  effect <- alpha[j] * alpha[k]

  effect_jacobian <- matrix(0, nrow = m, ncol = p)
  effect_jacobian[cbind(edge, j)] <- alpha[k]
  effect_jacobian[cbind(edge, k)] <- alpha[j]

  list(
    effect = effect,
    effect_jacobian = effect_jacobian
  )
}


.block_prepare <- function(
    mod,
    control_arr,
    diagnosed_arr,
    refresh_theta,
    est_n,
    nonpositive) {

  if (!requireNamespace("corrpops", quietly = TRUE)) {
    stop("The corrpops package is required.")
  }

  control_data <- corrpops::convert_corr_array_to_data_matrix(control_arr)
  diagnosed_data <- corrpops::convert_corr_array_to_data_matrix(diagnosed_arr)

  if (ncol(control_data) != ncol(diagnosed_data)) {
    stop("control_arr and diagnosed_arr must contain the same edges.")
  }

  m <- ncol(control_data)
  p <- as.integer(round((1 + sqrt(1 + 8 * m)) / 2))

  if (p * (p - 1) / 2 != m) {
    stop("The number of edges must equal p(p - 1)/2.")
  }

  alpha <- as.numeric(mod$alpha)
  theta_model <- as.numeric(mod$theta)

  if (length(alpha) != p || length(theta_model) != m) {
    stop("This implementation requires dim_alpha = 1 and vectorized theta.")
  }

  n_control <- nrow(control_data)
  n_diagnosed <- nrow(diagnosed_data)

  if (n_control < 2L || n_diagnosed < 2L) {
    stop("Each group must contain at least two subjects.")
  }

  pair_index <- which(lower.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
  effect_terms <- .block_effect_terms(alpha, pair_index)
  effect <- effect_terms$effect

  if (any(!is.finite(effect)) || any(effect == 0)) {
    stop("All multiplicative edge effects must be finite and nonzero.")
  }

  control_mean <- colMeans(control_data)
  diagnosed_mean <- colMeans(diagnosed_data)

  # The theta update used by the alternating corrpops estimator:
  #
  # theta(alpha) = [n_H rbar_H + n_D rbar_D / e(alpha)] / (n_H + n_D).
  theta_profile <- (
    n_control * control_mean +
      n_diagnosed * diagnosed_mean / effect
  ) / (n_control + n_diagnosed)

  theta <- if (refresh_theta) theta_profile else theta_model

  # Direct diagnosed-group derivative with theta held fixed:
  #
  # J_alpha = d[theta * e(alpha)] / d alpha
  #         = diag(theta) E,
  #
  # where E = d e(alpha) / d alpha.
  E <- effect_terms$effect_jacobian
  J_alpha <- sweep(E, 1L, theta, "*")
  mu_diagnosed <- theta * effect

  covariance_diagnosed <- corrpops::corrmat_covariance_from_datamatrix(
    diagnosed_data,
    est_n = est_n,
    nonpositive = nonpositive
  )
  W_diagnosed <- mod$weight_inv

  if (!is.null(mod$regularization$lambda) &&
      is.finite(mod$regularization$lambda) &&
      mod$regularization$lambda != 0) {
    warning(
      "The sandwich does not include derivatives of the alpha penalty.",
      call. = FALSE
    )
  }

  list(
    control_data = control_data,
    diagnosed_data = diagnosed_data,
    control_mean = control_mean,
    diagnosed_mean = diagnosed_mean,
    covariance_diagnosed = covariance_diagnosed,
    W_diagnosed = W_diagnosed,
    alpha = alpha,
    theta = theta,
    theta_model = theta_model,
    theta_profile = theta_profile,
    effect = effect,
    E = E,
    J_alpha = J_alpha,
    mu_diagnosed = mu_diagnosed,
    pair_index = pair_index,
    n_control = n_control,
    n_diagnosed = n_diagnosed,
    n_total = n_control + n_diagnosed,
    m = m,
    p = p
  )
}


.block_bread <- function(x) {
  # Joint estimating equations:
  #
  # U_theta = sum_H (r_i - theta)
  #         + sum_D (r_i / e(alpha) - theta)
  #
  # U_alpha = J_alpha' W_D sum_D [r_i - theta * e(alpha)].
  #
  # A = -d U / d(theta, alpha)'.

  # -d U_theta / d alpha'
  A_theta_alpha <- sweep(
    x$E,
    1L,
    colSums(x$diagnosed_data) / x$effect^2,
    "*"
  )

  diagnosed_residual_sum <-
    colSums(x$diagnosed_data) - x$n_diagnosed * x$mu_diagnosed
  weighted_residual_sum <-
    as.numeric(x$W_diagnosed %*% diagnosed_residual_sum)

  # -d U_alpha / d theta'
  A_alpha_theta <-
    x$n_diagnosed * crossprod(
      x$J_alpha,
      sweep(x$W_diagnosed, 2L, x$effect, "*")
    ) -
    t(sweep(x$E, 1L, weighted_residual_sum, "*"))

  # -d U_alpha / d alpha'.
  #
  # For edge (j,k), the only nonzero second derivatives of
  # mu_jk = theta_jk alpha_j alpha_k are
  #
  # d^2 mu_jk / (d alpha_j d alpha_k) = theta_jk.
  curvature <- matrix(0, nrow = x$p, ncol = x$p)
  j <- x$pair_index[, 1L]
  k <- x$pair_index[, 2L]
  edge_curvature <- x$theta * weighted_residual_sum
  curvature[cbind(j, k)] <- edge_curvature
  curvature[cbind(k, j)] <- edge_curvature

  A_alpha_alpha <-
    x$n_diagnosed * crossprod(
      x$J_alpha,
      x$W_diagnosed %*% x$J_alpha
    ) -
    curvature

  list(
    A_theta_theta = x$n_total,
    A_theta_alpha = A_theta_alpha,
    A_alpha_theta = A_alpha_theta,
    A_alpha_alpha = A_alpha_alpha
  )
}


.block_scores <- function(x, est_mu) {
  if (est_mu) {
    residual_control <- sweep(
      x$control_data,
      2L,
      x$theta,
      "-"
    )
    residual_diagnosed <- sweep(
      x$diagnosed_data,
      2L,
      x$mu_diagnosed,
      "-"
    )
  } else {
    residual_control <- sweep(
      x$control_data,
      2L,
      x$control_mean,
      "-"
    )
    residual_diagnosed <- sweep(
      x$diagnosed_data,
      2L,
      x$diagnosed_mean,
      "-"
    )
  }

  # One row per subject.
  theta_score_control <- residual_control
  theta_score_diagnosed <- sweep(
    residual_diagnosed,
    2L,
    x$effect,
    "/"
  )

  alpha_score_control <- matrix(
    0,
    nrow = x$n_control,
    ncol = x$p
  )
  alpha_score_diagnosed <-
    residual_diagnosed %*% x$W_diagnosed %*% x$J_alpha

  list(
    theta_control = theta_score_control,
    theta_diagnosed = theta_score_diagnosed,
    alpha_control = alpha_score_control,
    alpha_diagnosed = alpha_score_diagnosed
  )
}


.block_group_meat <- function(scores, finite_sample) {
  correction <- if (finite_sample) {
    nrow(scores) / (nrow(scores) - 1)
  } else {
    1
  }

  correction * crossprod(scores)
}


.block_symmetrize <- function(A) {
  (A + t(A)) / 2
}


.block_set_alpha_names <- function(A) {
  alpha_names <- paste0("alpha_", seq_len(nrow(A)))
  dimnames(A) <- list(alpha_names, alpha_names)
  A
}


# -----------------------------------------------------------------------------
# HC2 correction for the profiled p-dimensional scores
# -----------------------------------------------------------------------------

.block_symmetric_power <- function(A, power, eigenvalue_floor = 1e-10) {
  eig <- eigen(.block_symmetrize(A), symmetric = TRUE)
  values <- pmax(eig$values, eigenvalue_floor)^power
  eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
}


.block_hc2_multiplier <- function(total_bread, one_subject_bread) {
  B_half <- .block_symmetric_power(total_bread, 0.5)
  B_inverse_half <- .block_symmetric_power(total_bread, -0.5)

  leverage <- B_inverse_half %*%
    one_subject_bread %*%
    B_inverse_half
  leverage <- .block_symmetrize(leverage)

  eig <- eigen(leverage, symmetric = TRUE)
  lambda <- pmin(pmax(eig$values, 0), 1 - 1e-8)

  inverse_residual_half <- eig$vectors %*%
    diag((1 - lambda)^(-0.5), nrow = length(lambda)) %*%
    t(eig$vectors)

  B_inverse_half %*% inverse_residual_half %*% B_half
}


.block_apply_hc2 <- function(x, score_control, score_diagnosed,
                             est_n, nonpositive) {
  covariance_control <- corrpops::corrmat_covariance_from_datamatrix(
    x$control_data,
    est_n = est_n,
    nonpositive = nonpositive
  )
  W_control <- solve(covariance_control)

  # Derivatives of the two profiled means, theta(alpha) and
  # theta(alpha) * e(alpha).
  dtheta_dalpha <- -sweep(
    x$E,
    1L,
    (x$n_diagnosed / x$n_total) *
      x$diagnosed_mean / x$effect^2,
    "*"
  )

  D_control <- dtheta_dalpha
  D_diagnosed <- x$J_alpha + sweep(
    dtheta_dalpha,
    1L,
    x$effect,
    "*"
  )

  C_control <- crossprod(
    D_control,
    W_control %*% D_control
  )
  C_diagnosed <- crossprod(
    D_diagnosed,
    x$W_diagnosed %*% D_diagnosed
  )

  profile_bread <-
    x$n_control * C_control +
    x$n_diagnosed * C_diagnosed

  Q_control <- .block_hc2_multiplier(profile_bread, C_control)
  Q_diagnosed <- .block_hc2_multiplier(profile_bread, C_diagnosed)

  list(
    score_control = score_control %*% Q_control,
    score_diagnosed = score_diagnosed %*% Q_diagnosed,
    covariance_control = covariance_control,
    profile_bread = profile_bread,
    C_control = C_control,
    C_diagnosed = C_diagnosed
  )
}


# -----------------------------------------------------------------------------
# Public functions
# -----------------------------------------------------------------------------

#' Compute a Dense Joint-Block GEE Sandwich Variance
#'
#' Estimates the covariance matrix of the fitted alpha parameters using a
#' Generalized Estimating Equations (GEE) sandwich estimator for the joint
#' parameter vector \eqn{(\theta, \alpha)}.
#'
#' This function forms the complete joint sandwich estimator for the theta and
#' alpha estimating equations and returns the alpha-alpha covariance block. In
#' contrast to an alpha-only sandwich, the calculation accounts for the
#' variability introduced by estimating theta.
#'
#' The function is intended primarily as a transparent reference
#' implementation. If \eqn{m = p(p - 1)/2} is the number of correlation edges,
#' it explicitly forms and inverts a bread matrix of dimension
#' \eqn{(m + p) \times (m + p)}. For larger problems,
#' \code{\link{compute_schur_block_variance}} gives the same alpha covariance
#' without forming the full joint matrix.
#'
#' @param mod A fitted model object returned by
#'   \code{\link[corrpops]{estimate_model}}. The implementation requires a
#'   one-dimensional multiplicative-identity model with a length-p alpha
#'   vector and a length-m vectorized theta parameter.
#' @param control_arr Correlation matrices from the control group. This may be
#'   an array of correlation matrices or a matrix in which each row contains a
#'   vectorized correlation matrix, as accepted by
#'   \code{\link[corrpops]{convert_corr_array_to_data_matrix}}.
#' @param diagnosed_arr Correlation matrices from the diagnosed group, in the
#'   same format and with the same correlation edges as \code{control_arr}.
#' @param est_mu Logical. If \code{TRUE}, subject-level scores are computed from
#'   residuals around the fitted group means. If \code{FALSE}, residuals are
#'   centered around the corresponding sample group means. Default is
#'   \code{TRUE}.
#' @param refresh_theta Logical. If \code{TRUE}, theta is recomputed from the
#'   fitted alpha using the theta update of the alternating corrpops estimator.
#'   If \code{FALSE}, the value stored in \code{mod$theta} is used. Default is
#'   \code{TRUE}.
#' @param est_n Logical argument passed to
#'   \code{corrpops::corrmat_covariance_from_datamatrix()} when estimating the
#'   diagnosed-group working covariance. Default is \code{TRUE}.
#' @param nonpositive Character argument passed to
#'   \code{corrpops::corrmat_covariance_from_datamatrix()} to determine how a
#'   non-positive-definite working covariance is handled. Default is
#'   \code{"force"}.
#' @param finite_sample Logical. If \code{TRUE}, each group's empirical meat
#'   contribution is multiplied by \eqn{n_g/(n_g - 1)}. Default is
#'   \code{TRUE}.
#' @param return_details Logical. If \code{FALSE}, return only the estimated
#'   covariance matrix of alpha. If \code{TRUE}, return intermediate matrices
#'   and fitted quantities in addition to the covariance matrix. Default is
#'   \code{FALSE}.
#'
#' @details
#' Let \eqn{U(\theta, \alpha)} denote the stacked theta and alpha estimating
#' equations. The estimator has the sandwich form
#' \deqn{\widehat{\operatorname{Var}}(\widehat\eta)
#' = A^{-1} B A^{-T},}
#' where \eqn{\eta = (\theta, \alpha)}, \eqn{A} is the negative derivative
#' of the stacked estimating equations, and \eqn{B} is the empirical
#' cross-product of the subject-level estimating functions. The returned
#' covariance matrix is the alpha-alpha block of this joint sandwich.
#'
#' Derivatives of an alpha regularization penalty are not included. A warning
#' is issued when the fitted model contains a nonzero regularization parameter.
#'
#' @return
#' If \code{return_details = FALSE}, a symmetric \eqn{p \times p} matrix
#' containing the estimated covariance matrix of \code{mod$alpha}.
#'
#' If \code{return_details = TRUE}, a list with components:
#' \describe{
#'   \item{variance}{The alpha-alpha covariance block.}
#'   \item{joint_variance}{The complete covariance matrix for
#'     \eqn{(\theta, \alpha)}.}
#'   \item{bread}{The complete joint bread matrix.}
#'   \item{meat}{The complete joint empirical meat matrix.}
#'   \item{score_control, score_diagnosed}{Matrices of subject-level joint
#'     scores for the two groups.}
#'   \item{blocks}{The four blocks of the joint bread matrix.}
#'   \item{alpha}{The fitted alpha vector.}
#'   \item{theta}{The theta vector used in the calculation.}
#'   \item{theta_model}{The theta vector stored in \code{mod}.}
#'   \item{theta_profile}{Theta recomputed from the fitted alpha and the two
#'     group sample means.}
#'   \item{working_covariance_diagnosed}{The estimated diagnosed-group working
#'     covariance matrix.}
#' }
#'
#' @seealso
#' \code{\link{compute_schur_block_variance}},
#' \code{\link[corrpops]{estimate_model}},
#' \code{\link[corrpops]{convert_corr_array_to_data_matrix}}
#'
#' @export
#'
compute_dense_block_variance <- function(
    mod,
    control_arr,
    diagnosed_arr,
    est_mu = TRUE,
    refresh_theta = TRUE,
    est_n = TRUE,
    nonpositive = "force",
    finite_sample = TRUE,
    return_details = FALSE) {

  x <- .block_prepare(
    mod = mod,
    control_arr = control_arr,
    diagnosed_arr = diagnosed_arr,
    refresh_theta = refresh_theta,
    est_n = est_n,
    nonpositive = nonpositive
  )
  A <- .block_bread(x)
  scores <- .block_scores(x, est_mu = est_mu)

  bread <- rbind(
    cbind(x$n_total * diag(x$m), A$A_theta_alpha),
    cbind(A$A_alpha_theta, A$A_alpha_alpha)
  )

  score_control <- cbind(
    scores$theta_control,
    scores$alpha_control
  )
  score_diagnosed <- cbind(
    scores$theta_diagnosed,
    scores$alpha_diagnosed
  )

  meat <-
    .block_group_meat(score_control, finite_sample) +
    .block_group_meat(score_diagnosed, finite_sample)

  bread_inverse <- solve(bread)
  joint_variance <- bread_inverse %*% meat %*% t(bread_inverse)
  joint_variance <- .block_symmetrize(joint_variance)

  alpha_index <- x$m + seq_len(x$p)
  variance <- joint_variance[alpha_index, alpha_index, drop = FALSE]
  variance <- .block_set_alpha_names(variance)

  if (!return_details) {
    return(variance)
  }

  list(
    variance = variance,
    joint_variance = joint_variance,
    bread = bread,
    meat = meat,
    score_control = score_control,
    score_diagnosed = score_diagnosed,
    blocks = A,
    alpha = x$alpha,
    theta = x$theta,
    theta_model = x$theta_model,
    theta_profile = x$theta_profile,
    working_covariance_diagnosed = x$covariance_diagnosed
  )
}


#' Compute a GEE Sandwich Variance via Structured Schur-Block GEE Sandwich Variance
#'
#' Estimates the covariance matrix of the fitted alpha parameters using the
#' same joint theta-alpha estimating equations as
#' \code{\link{compute_dense_block_variance}}, while eliminating the theta
#' block through a Schur complement.
#'
#' This is the computationally efficient implementation of the joint-block GEE
#' sandwich estimator. It accounts for variability in theta estimation but
#' avoids forming or inverting the full \eqn{(m + p) \times (m + p)} bread
#' matrix.
#'
#' @inheritParams compute_dense_block_variance
#' @param hc2 Logical. If \code{TRUE}, apply an HC2 leverage correction to the
#'   profiled, p-dimensional subject scores before constructing the meat.
#'   Default is \code{FALSE}.
#'
#' @details
#' The theta-theta bread block is \eqn{(n_H + n_D) I_m}. The function uses this
#' structure to form the Schur-complement bread
#' \deqn{A_{\alpha\alpha \cdot \theta}
#' = A_{\alpha\alpha}
#' - A_{\alpha\theta} A_{\theta\theta}^{-1}
#'   A_{\theta\alpha}}
#' and the corresponding efficient subject scores
#' \deqn{u_{\alpha i}^{\mathrm{eff}}
#' = u_{\alpha i}
#' - A_{\alpha\theta} A_{\theta\theta}^{-1} u_{\theta i}.}
#' The resulting sandwich requires inversion only of a \eqn{p \times p}
#' matrix and, without HC2, is algebraically equivalent to extracting the
#' alpha-alpha block from the dense joint sandwich.
#'
#' When \code{hc2 = TRUE}, the function computes group-specific leverage
#' adjustments using the profiled alpha bread and applies the resulting matrix
#' multipliers to the effective subject scores.
#'
#' Derivatives of an alpha regularization penalty are not included. A warning
#' is issued when the fitted model contains a nonzero regularization parameter.
#'
#' @return
#' If \code{return_details = FALSE}, a symmetric \eqn{p \times p} matrix
#' containing the estimated covariance matrix of \code{mod$alpha}.
#'
#' If \code{return_details = TRUE}, a list with components:
#' \describe{
#'   \item{variance}{The estimated covariance matrix of alpha.}
#'   \item{schur_bread}{The p-dimensional Schur-complement bread matrix.}
#'   \item{meat}{The empirical meat formed from the effective alpha scores.}
#'   \item{effective_score_control, effective_score_diagnosed}{Matrices of
#'     subject-level alpha scores after eliminating theta.}
#'   \item{alpha_from_theta}{The matrix
#'     \eqn{A_{\alpha\theta} A_{\theta\theta}^{-1}}.}
#'   \item{blocks}{The four blocks of the original joint bread matrix.}
#'   \item{alpha}{The fitted alpha vector.}
#'   \item{theta}{The theta vector used in the calculation.}
#'   \item{theta_model}{The theta vector stored in \code{mod}.}
#'   \item{theta_profile}{Theta recomputed from the fitted alpha and the two
#'     group sample means.}
#'   \item{working_covariance_diagnosed}{The estimated diagnosed-group working
#'     covariance matrix.}
#'   \item{hc2}{Whether the HC2 correction was requested.}
#'   \item{hc2_details}{If HC2 was requested, its profile bread, group bread
#'     contributions, control working covariance, and adjusted scores;
#'     otherwise \code{NULL}.}
#' }
#'
#' @seealso
#' \code{\link{compute_dense_block_variance}},
#' \code{\link[corrpops]{estimate_model}},
#' \code{\link[corrpops]{convert_corr_array_to_data_matrix}}
#'
#' @export
#'
compute_schur_variance <- function(
    mod,
    control_arr,
    diagnosed_arr,
    est_mu = TRUE,
    hc2 = FALSE,
    refresh_theta = TRUE,
    est_n = TRUE,
    nonpositive = "force",
    finite_sample = TRUE,
    return_details = FALSE) {

  x <- .block_prepare(
    mod = mod,
    control_arr = control_arr,
    diagnosed_arr = diagnosed_arr,
    refresh_theta = refresh_theta,
    est_n = est_n,
    nonpositive = nonpositive
  )
  A <- .block_bread(x)
  scores <- .block_scores(x, est_mu = est_mu)

  # A_alpha_theta A_theta_theta^{-1}; the theta bread is n_total I_m.
  alpha_from_theta <- A$A_alpha_theta / x$n_total

  schur_bread <-
    A$A_alpha_alpha -
    alpha_from_theta %*% A$A_theta_alpha

  # Efficient scores after eliminating theta:
  #
  # u_alpha^eff = u_alpha -
  #               A_alpha_theta A_theta_theta^{-1} u_theta.
  effective_score_control <-
    scores$alpha_control -
    scores$theta_control %*% t(alpha_from_theta)

  effective_score_diagnosed <-
    scores$alpha_diagnosed -
    scores$theta_diagnosed %*% t(alpha_from_theta)

  hc2_details <- NULL
  if (hc2) {
    hc2_details <- .block_apply_hc2(
      x = x,
      score_control = effective_score_control,
      score_diagnosed = effective_score_diagnosed,
      est_n = est_n,
      nonpositive = nonpositive
    )
    effective_score_control <- hc2_details$score_control
    effective_score_diagnosed <- hc2_details$score_diagnosed
  }

  meat <-
    .block_group_meat(effective_score_control, finite_sample) +
    .block_group_meat(effective_score_diagnosed, finite_sample)

  schur_inverse <- solve(schur_bread)
  variance <- schur_inverse %*% meat %*% t(schur_inverse)
  variance <- .block_set_alpha_names(.block_symmetrize(variance))

  if (!return_details) {
    return(variance)
  }

  list(
    variance = variance,
    schur_bread = schur_bread,
    meat = meat,
    effective_score_control = effective_score_control,
    effective_score_diagnosed = effective_score_diagnosed,
    alpha_from_theta = alpha_from_theta,
    blocks = A,
    alpha = x$alpha,
    theta = x$theta,
    theta_model = x$theta_model,
    theta_profile = x$theta_profile,
    working_covariance_diagnosed = x$covariance_diagnosed,
    hc2 = hc2,
    hc2_details = hc2_details
  )
}


# Compatibility aliases for earlier simulation scripts.

#' @rdname compute_schur_block_variance
#' @export
compute_structured_schur_variance <- compute_schur_variance

#' @rdname compute_dense_block_variance
#' @export
compute_full_block_variance <- compute_dense_block_variance
