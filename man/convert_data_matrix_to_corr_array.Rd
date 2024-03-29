% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_corr_array_to_data_matrix.R
\name{convert_data_matrix_to_corr_array_raw}
\alias{convert_data_matrix_to_corr_array_raw}
\alias{convert_data_matrix_to_corr_array}
\title{Convert Data Matrix of Correlation Matrices to an Array}
\usage{
convert_data_matrix_to_corr_array_raw(matr)

convert_data_matrix_to_corr_array(matr, verbose = FALSE)
}
\arguments{
\item{matr}{a data matrix of vectorized correlation matrices}

\item{verbose}{should add messages?}
}
\value{
an array with the rows as correlation matrices
}
\description{
Convert a data matrix of vectorized correlation matrices to an array.
Each row becomes a matrix in the output array. If matr is already of class array, nothing is performed.
}
\section{Functions}{
\itemize{
\item \code{convert_data_matrix_to_corr_array_raw}: internal function of convert_data_matrix_to_corr_array
}}

\seealso{
Other vectriangle: 
\code{\link{convert_corr_array_to_data_matrix_raw}()},
\code{\link{triangle2vector}()}
}
\concept{vectriangle}
