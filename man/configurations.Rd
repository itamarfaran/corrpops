% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_configs.R
\name{configurations}
\alias{configurations}
\title{Default Configurations for Model Estimation}
\usage{
configurations(index)
}
\arguments{
\item{index}{which configurations to print. if missing, print all. can be one or more of the aforementioned.}
}
\description{
A list of default configurations for \link[corrpops]{estimate_model}. running the function will print the defaults.
\itemize{
\item model_reg: a list with params lambda, lp for regularization of alpha. the loss function has an element of reg_lambda * sum(abs(alpha - LinkFunc$null_value)^reg_p) ^ (1/reg_p)
\item matrix_reg: parameters passed on to \link[corrpops]{regularize_matrix}. the weighting matrix (covariance matrix of correlations) is regularized accordingly
\item iterations: parameters used for the iterations' stopping rule. the definitions of 'maxit', 'reltol', 'abstol' are the same is in \link[stats]{optim}. iterations will stop only if there were 'minit' optim iterations with convergence=0 in a row.
\item optim: possible parameters to be passed to \link[stats]{optim}, namely 'method', 'reltol', 'abstol'. also, log_optim is a boolean - whether to save the result of the call to \link[stats]{optim} in each iteration.
}
}
