% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrmat_covariance.R
\name{corrmat_covariance}
\alias{corrmat_covariance}
\title{Calculate the Covariance Matrix of a Correlation Matrix}
\usage{
corrmat_covariance(
  matr,
  fisher_z = FALSE,
  nonpositive = c("stop", "force", "ignore"),
  use_cpp = TRUE
)
}
\arguments{
\item{matr}{the correlation matrix. can be vectorized from \link[corrpops]{triangle2vector}}

\item{fisher_z}{if true, calculate the covariance matrix of a fisher-Z transformed correlation matrix. It is assumed that the correlations are already Fisher transformed, and the matrix will under go the Inverse Fisher Transformation.}

\item{nonpositive}{error handling when matrix is not positive definite. can be one of 'stop', 'force', 'ignore'. if 'force' is chosen, \link[Matrix]{nearPD} will be used.}

\item{use_cpp}{whether to use c++ source code. default true.}
}
\value{
the covariance matrix
}
\description{
Calculate the covariance matrix of a correlation matrix. The matrix can be in vectorized form.
}
\seealso{
Other corrcalc: 
\code{\link{corrcalc_c}()},
\code{\link{corrmat_covariance_from_datamatrix}()}
}
\concept{corrcalc}
