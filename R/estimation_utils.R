theta_of_alpha <- function(alpha, control_datamatrix, diagnosed_datamatrix, LinkFunc, d = 1){
  out <- rbind(LinkFunc$rev_func(datamatrix = diagnosed_datamatrix, a = alpha, d = d), control_datamatrix)
  out <- colMeans(out)
  return(out)
}


sum_of_squares <- function(
  alpha, theta, diagnosed_datamatrix, inv_sigma, LinkFunc,
  sigma, dim_alpha = 1, reg_lambda = 0, reg_p = 2)
{
  if(missing(inv_sigma))
    inv_sigma <- solve(sigma)

  g11 <- as.matrix(triangle2vector(LinkFunc$func(t = theta, a = alpha, d = dim_alpha)))
  sse <- nrow(diagnosed_datamatrix) * t(g11) %*% inv_sigma %*% ( 0.5 * g11 - colMeans(diagnosed_datamatrix) )

  if(reg_lambda > 0)
    sse <- sse + reg_lambda*sum((alpha - LinkFunc$null_value)^reg_p)

  return(sse)
}


get_update_message <- function(i, start_time, convergence, distance){
  msg <- paste0(i, " (", round(as.double.difftime(Sys.time() - start_time, units = "secs")),
                "s, ", convergence, ", ", round(distance, 5), "); ")
  return(msg)
}


optimiser <- function(
  control_datamatrix, diagnosed_datamatrix, alpha0, theta0, weight_matrix, dim_alpha, LinkFunc,
  model_reg_config, matrix_reg_config, iter_config, optim_config, early_stop, verbose)
{

  if('reltol' %in% names(iter_config) & 'abstol' %in% names(iter_config))
    stop('can supply only one of reltol or abstol')

  model_reg_config <- utils::modifyList(default_model_reg_config, model_reg_config)
  matrix_reg_config <- utils::modifyList(default_matrix_reg_config, matrix_reg_config)
  iter_config <- utils::modifyList(default_iter_config, iter_config)
  optim_config <- utils::modifyList(default_optim_config, optim_config)

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

  for(i in 2:iter_config$max_loop){
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
    if(i > iter_config$min_loop){
      look_back <- iter_config$min_loop - 1
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

  if(i == iter_config$max_loop)
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
    vcov = weight_matrix_reg_inv,
    convergence = convergence,
    steps = steps, log_optim = log_optim_out
  )

  return(output)
}
