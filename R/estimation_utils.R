#' Compute Theta as a Function of Alpha
#'
#' @param alpha estimate of alpha
#' @param control_datamatrix data matrix of control subjects
#' @param diagnosed_datamatrix data matrix of diagnosed subjects
#' @param LinkFunc \link[corrpops]{LinkFuncSkeleton}
#' @param d dimension of alpha
#' @return estimate of theta based on control and doagnosed subjects
theta_of_alpha <- function(alpha, control_datamatrix, diagnosed_datamatrix, LinkFunc, d = 1){
  reversed_diagnosed_datamatrix <- LinkFunc$rev_func(datamatrix = diagnosed_datamatrix, a = alpha, d = d)
  joint_datamatrix <- rbind(reversed_diagnosed_datamatrix, control_datamatrix)
  out <- colMeans(joint_datamatrix)
  return(out)
}


#' Loss Function for Optimization Process on Alpha
#'
#' @param alpha estimate of alpha
#' @param theta estimate of theta
#' @param diagnosed_datamatrix data matrix of diagnosed subjects
#' @param inv_sigma covariance matrix of correlation data matrix, inverted
#' @param LinkFunc \link[corrpops]{LinkFuncSkeleton}
#' @param sigma if inv_sigma is missing, invert this matrix and use it as inv_sigma
#' @param dim_alpha dimension of alpha
#' @param reg_lambda see model_reg in \link[corrpops]{configurations}
#' @param reg_p see model_reg in \link[corrpops]{configurations}
#' @return value of sum of squares
sum_of_squares <- function(
  alpha, theta, diagnosed_datamatrix, inv_sigma, LinkFunc,
  sigma, dim_alpha = 1, reg_lambda = 0, reg_p = 2)
{
  if(missing(inv_sigma))
    inv_sigma <- solve(sigma)

  g11 <- as.matrix(triangle2vector(LinkFunc$func(t = theta, a = alpha, d = dim_alpha)))
  sse <- nrow(diagnosed_datamatrix) * t(g11) %*% inv_sigma %*% ( 0.5 * g11 - colMeans(diagnosed_datamatrix) )

  if(reg_lambda > 0)
    sse <- sse + reg_lambda * sum(abs(alpha - LinkFunc$null_value)^reg_p) ^ (1/reg_p)

  return(sse)
}


get_update_message <- function(i, start_time, convergence, distance){
  msg <- paste0(i, " (", round(as.double.difftime(Sys.time() - start_time, units = "secs")),
                "s, ", convergence, ", ", round(distance, 5), "); ")
  return(msg)
}


#' Iterative Call for optim
#'
#' @param control_datamatrix data matrix of control subjects
#' @param diagnosed_datamatrix data matrix of diagnosed subjects
#' @param alpha0 starting point for alpha
#' @param theta0 starting point for theta
#' @param weight_matrix the weighting matrix to be used in the mahalonobis sum of squares
#' @param dim_alpha dimension of alpha
#' @param LinkFunc \link[corrpops]{LinkFuncSkeleton}
#' @param model_reg_config see \link[corrpops]{configurations}
#' @param matrix_reg_config see \link[corrpops]{configurations}
#' @param iter_config see \link[corrpops]{configurations}
#' @param optim_config see \link[corrpops]{configurations}
#' @param early_stop if true, stop the optimization of the joint loss function (of theta and alpha) didn't decrease in the last iteration.
#' @param verbose if true, print status to console
#' @return see \link[corrpops]{estimate_model}
#' @seealso \link[corrpops]{estimate_model}
optimiser <- function(
  control_datamatrix, diagnosed_datamatrix, alpha0, theta0, weight_matrix, dim_alpha, LinkFunc,
  model_reg_config, matrix_reg_config, iter_config, optim_config, early_stop, verbose)
{

  if('reltol' %in% names(iter_config) & 'abstol' %in% names(iter_config))
    stop('can supply only one of reltol or abstol')

  model_reg_config <- utils::modifyList(configs$model_reg, model_reg_config)
  matrix_reg_config <- utils::modifyList(configs$matrix_reg, matrix_reg_config)
  iter_config <- utils::modifyList(configs$iterations, iter_config)
  optim_config <- utils::modifyList(configs$optim, optim_config)

  p <- .5 * (1 + sqrt(1 + 8 * ncol(diagnosed_datamatrix)))
  m <- .5 * p * (p - 1)

  if(is.null(theta0))
    theta0 <- colMeans(rbind(control_datamatrix, diagnosed_datamatrix))
  if(is.null(alpha0))
    alpha0 <- matrix(LinkFunc$null_value, nrow = p, ncol = dim_alpha)

  dim_alpha <- length(alpha0) / p
  if(dim_alpha %% 1 != 0)
    stop("alpha0 not multiplicative of p")

  if(!matrixcalc::is.positive.semi.definite(vector2triangle(theta0, diag_value = 1)) ||
      !matrixcalc::is.positive.semi.definite(LinkFunc$func(t = theta0, a = alpha0, d = dim_alpha)))
    warning("Initial parameters dont result with positive-definite matrices")

  if(is.null(weight_matrix)) {
    weight_matrix <- weight_matrix_reg <- weight_matrix_reg_inv <- diag(m)
  } else {
    weight_matrix_reg <- regularize_matrix(
      weight_matrix,
      method = matrix_reg_config$method,
      const = matrix_reg_config$const
    )
    weight_matrix_reg_inv <- solve(weight_matrix_reg)  # todo: bottleneck
  }

  temp_theta <- theta0
  temp_alpha <- alpha0
  log_optim_out <- list()
  steps <- list()
  convergence <- NA_integer_

  steps[[1]] <- list(
    theta = temp_theta,
    alpha = temp_alpha,
    value = sum_of_squares(
      theta = temp_theta,
      alpha = temp_alpha,
      diagnosed_datamatrix = diagnosed_datamatrix,
      inv_sigma = weight_matrix_reg_inv,
      LinkFunc = LinkFunc,
      dim_alpha = dim_alpha,
      reg_lambda = model_reg_config$lambda,
      reg_p = model_reg_config$lp,
    )
  )

  start_time <- Sys.time()
  if(verbose)
    message(paste0("Time of intialization: ", start_time, "; Progress: 'Loop, (Time, Convergence, Distance)'"))

  for(i in 2:iter_config$maxit){
    temp_theta <- theta_of_alpha(
      alpha = temp_alpha,
      control_datamatrix = control_datamatrix,
      diagnosed_datamatrix = diagnosed_datamatrix,
      LinkFunc = LinkFunc,
      d = dim_alpha)

    optim_alpha <- stats::optim(
      par = temp_alpha,
      fn = sum_of_squares,
      theta = temp_theta,
      diagnosed_datamatrix = diagnosed_datamatrix,
      inv_sigma = weight_matrix_reg_inv,
      LinkFunc = LinkFunc,
      dim_alpha = dim_alpha,
      reg_lambda = model_reg_config$lambda,
      reg_p = model_reg_config$lp,
      method = optim_config$method,
      control = list(
        maxit = min(max(500, i*100), 2000),
        reltol = optim_config$reltol
      )
    )

    temp_alpha <- optim_alpha$par
    steps[[i]] <- list(
      theta = temp_theta,
      alpha = temp_alpha,
      value = optim_alpha$value,
      convergence = optim_alpha$convergence
    )
    convergence <- c(convergence, optim_alpha$convergence)
    log_optim_out[[i]] <- if(optim_config$log_optim) optim_alpha else NA

    # Stopping rule
    if('abstol' %in% names(iter_config)){
      distance <- vect_norm(steps[[i]]$alpha - steps[[i-1]]$alpha, sqrt = FALSE)
      distance_lower_than_threshold <-
        distance < iter_config$abstol
    } else {
      distance <- abs(steps[[i - 1]]$value - steps[[i]]$value)
      distance_lower_than_threshold <-
        distance < (iter_config$reltol * (abs(steps[[i]]$value) + iter_config$reltol))
    }

    if(verbose)
      cat(get_update_message(i, start_time, steps[[i]]$convergence, distance))

    stopping_condition <- FALSE
    if(i > iter_config$minit){
      look_back <- iter_config$minit - 1
      index <- if(look_back > 0) i - 0:look_back else i
      stopping_condition <- distance_lower_than_threshold & (sum(convergence[index]) == 0)
    }

    if(early_stop){
      did_converge <- convergence[length(convergence)] == 0
      did_minimize <- steps[[i]]$value <= steps[[i-1]]$value

      if(did_converge & !did_minimize){
        steps[[i]] <- NULL
        i <- i - 1

        temp_theta <- steps[[i]]$theta
        temp_alpha <- steps[[i]]$alpha
        convergence[length(convergence)] <- -1
        stopping_condition <- TRUE
        warning('early stopping used; last iteration didn\'t minimize target')
      }
    }

    if(stopping_condition)
      break()
  }

  if(i == iter_config$maxit)
    warning('optimization reached maximum iterations')

  if(verbose){
    total_time <- Sys.time() - start_time
    units(total_time) <- 'secs'
    total_time <- as.numeric(total_time)
    message(paste0("\nTotal time: ", floor(total_time/60), " minutes and ", round(total_time %% 60, 1), " seconds."))
  }

  output <- list(
    theta = temp_theta,
    alpha = temp_alpha,
    LinkFunc = LinkFunc,
    regularization = model_reg_config,
    vcov = weight_matrix_reg_inv,
    convergence = convergence,
    steps = steps,
    log_optim = log_optim_out
  )

  return(output)
}
