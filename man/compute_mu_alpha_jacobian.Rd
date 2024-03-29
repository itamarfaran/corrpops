% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_gee_variance.R
\name{compute_mu_alpha_jacobian}
\alias{compute_mu_alpha_jacobian}
\title{Compute the Jacobian of the link function by Alpha parameter}
\usage{
compute_mu_alpha_jacobian(
  group,
  alpha,
  control_datamatrix,
  diagnosed_datamatrix,
  d = 1,
  LinkFunc
)
}
\arguments{
\item{group}{can be either 'diagnosed' or 'control'}

\item{alpha}{the vector to calculate the Jacobian on}

\item{control_datamatrix}{array of correlation matrices of control group, vectorized and organized in a data matrix}

\item{diagnosed_datamatrix}{same as control_datamatrix but of the diagnosed group}

\item{d}{the number of columns in alpha, default 1}

\item{LinkFunc}{the link function $g$ to use}
}
\value{
the Jacobian matrix
}
\description{
Helper function to calculate the Jacobian of the expected value by Alpha:
\eqn{J_{\alpha}\left(g\left(\Theta,\alpha\right)\right)_{ij}=
\frac{\partial g_{i}\left(\Theta,\alpha\right)}{\partial\alpha_{j}}}
}
