theta_of_alpha <- function(alpha, healthy_dt, sick_dt, linkFun, d = 1){
  out <- rbind(linkFun$CLEAN(dt = sick_dt, a = alpha, d = d), healthy_dt)
  out <- colMeans(out)
  return(out)
}


sum_of_squares <- function(
  alpha, theta, sick_dt, inv_sigma,
  linkFun, sigma, dim_alpha = 1,
  reg_lambda = 0, reg_p = 2)
  {

  if(missing(inv_sigma))
    inv_sigma <- solve(sigma)

  g11 <- as.matrix(triangle2vector(linkFun$FUN(t = theta, a = alpha, d = dim_alpha)))
  sse <- nrow(sick_dt) * t(g11) %*% inv_sigma %*% ( 0.5 * g11 - colMeans(sick_dt) )

  if(reg_lambda > 0)
    sse <- sse + reg_lambda*sum((alpha - linkFun$NULL_VAL)^reg_p)
  return(sse)
}


inner_optim_loop <- function(
  healthy_dt, sick_dt, alpha0 = NULL, theta0 = NULL, weight_matrix = NULL,
  dim_alpha = 1, linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iter_config = list(), optim_config = list(), early_stop = FALSE,
  verbose = TRUE){

  if('reltol' %in% names(iter_config) & 'abstol' %in% names(iter_config))
    stop('can supply only one of reltol or abstol')

  model_reg_config <- modifyList(list(lambda = 0, lp = 2), model_reg_config)
  matrix_reg_config <- modifyList(list(method = 'constant', const = 0), matrix_reg_config)
  iter_config <- modifyList(list(max_loop = 50, reltol = 1e-06, min_loop = 3), iter_config)
  optim_config <- modifyList(list(method = "BFGS", reltol = 1e-06, log_optim = FALSE), optim_config)

  p <- 0.5 + sqrt(1 + 8*ncol(sick_dt))/2
  m <- 0.5*p*(p-1)

  # if(is.null(theta0)) theta0 <- colMeans(healthy_dt)
  if(is.null(theta0)) theta0 <- colMeans(rbind(healthy_dt, sick_dt))
  if(is.null(alpha0)) alpha0 <- matrix(linkFun$NULL_VAL, nr = p, nc = dim_alpha)

  if(!is.positive.semi.definite(vector2triangle(theta0, diag_value = 1))
     || !is.positive.semi.definite(linkFun$FUN(t = theta0, a = alpha0, d = dim_alpha)))
    warning("Initial parameters dont result with positive-definite matrices")

  dim_alpha <- length(alpha0)/p
  if(dim_alpha %% 1 != 0) stop("alpha0 not multiplicative of p")

  if(is.null(weight_matrix)) {
    weight_matrix <- weight_matrix_reg <- weight_matrix_reg_inv <- diag(m)
  } else {
    weight_matrix_reg <- regularize_matrix(
      weight_matrix,
      method = matrix_reg_config$method,
      const = matrix_reg_config$const
    )
    weight_matrix_reg_inv <- solve(weight_matrix_reg)
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
      sick_dt = sick_dt,
      inv_sigma = weight_matrix_reg_inv,
      linkFun = linkFun,
      dim_alpha = dim_alpha,
      reg_lambda = model_reg_config$lambda,
      reg_p = model_reg_config$lp,
    )
  )


  tt <- Sys.time()
  if(verbose) message(paste0("Time of intialization: ", tt, "; Progress: 'Loop, (Time, Convergence, Distance)'"))
  for(i in 2:iter_config$max_loop){

    temp_theta <- theta_of_alpha(
      alpha = temp_alpha,
      healthy_dt = healthy_dt,
      sick_dt = sick_dt,
      linkFun = linkFun,
      d = dim_alpha)

    optim_alpha <- optim(
      par = temp_alpha,
      fn = sum_of_squares,
      theta = temp_theta,
      sick_dt = sick_dt,
      inv_sigma = weight_matrix_reg_inv,
      linkFun = linkFun,
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
      distance <- sqrt(mean((steps[[i]]$alpha - steps[[i-1]]$alpha)^2))
      distance_lower_than_threshold <-  distance < iter_config$abstol
    } else {
      distance <- abs(steps[[i-1]]$value - steps[[i]]$value)
      distance_lower_than_threshold <- distance < (iter_config$reltol * (abs(steps[[i]]$value) + iter_config$reltol))
    }

    if(verbose) cat(paste0(
      i, " (", round(as.double.difftime(Sys.time() - tt, units = "secs")), "s, ",
      steps[[i]]$convergence, ", ", round(distance, 5), "); "
    ))

    condition0 <- FALSE
    if(i > iter_config$min_loop){
      look_back <- iter_config$min_loop - 1
      index <- if(look_back > 0) i - 0:look_back else i
      condition0 <- distance_lower_than_threshold & (sum(convergence[index]) == 0)
    }

    # Early Stop
    if(early_stop){
      did_converge <- convergence[length(convergence)] == 0
      did_minimize <- steps[[i]]$value <= steps[[i-1]]$value

      if(did_converge & !did_minimize){
        steps[[i]] <- NULL
        i <- i - 1

        temp_theta <- steps[[i]]$theta
        temp_alpha <- steps[[i]]$alpha
        convergence[length(convergence)] <- -1
        condition0 <- TRUE
        warning('early stopping used; last iteration didn\'t minimize target')
      }
    }

    if(condition0) break()
  }
  if(i == iter_config$max_loop) warning('optimization reached maximum iterations')
  if(verbose){
    tt <- Sys.time() - tt
    units(tt) <- "secs"
    tt <- as.numeric(tt)
    message(paste0("\nTotal time: ", floor(tt/60), " minutes and ", round(tt %% 60, 1), " seconds."))
  }

  output <- list(
    theta = temp_theta,
    alpha = temp_alpha,
    linkFun = linkFun,
    vcov = weight_matrix_reg_inv,
    convergence = convergence,
    steps = steps, log_optim = log_optim_out
  )

  return(output)
}
