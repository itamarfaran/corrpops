% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_jackknife.R
\name{infer_jackknife}
\alias{infer_jackknife}
\title{Infer Jackknife Results}
\usage{
infer_jackknife(obj)
}
\arguments{
\item{obj}{output object of \link[corrpops]{estimate_model_jacknife}}
}
\value{
a list consisting of estimates and variance matrix

\itemize{
\item estimate: Jacknife estimate of alpha (average estimate over all jackknife iterations)
}

\itemize{
\item variance: Jacknife estimate of the variance matrix
}
}
\description{
Calculate the Jackknife Variance based on Jackknife Estimation
}
