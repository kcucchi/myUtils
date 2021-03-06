% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_unique_cubic_spline.R
\name{dist_param}
\alias{dist_param}
\title{Definition of the distance function used in optimisation}
\usage{
dist_param(param, t, y)
}
\arguments{
\item{param}{the vector of parameters to be optimized (amplitude, timing and width)}

\item{t}{numeric vector of time values}

\item{y}{numeric vector of observed values at time \code{t}}
}
\value{
the residual sum of squares between N_4(t,x_0,A,sigma) and y
}
\description{
Definition of the distance function used in optimisation
}
