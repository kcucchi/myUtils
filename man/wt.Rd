% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wavelet_morlet.R
\name{wt}
\alias{wt}
\title{Perform Morlet Wavelet transform}
\usage{
wt(y, J = 110, vect_time = NULL, dt = 1)
}
\arguments{
\item{y}{timeseries to transform}

\item{J}{determines the largest scale s_J = s_0 2^(J dj)}

\item{vect_time}{(optional) vector ot time points corresponding to y}

\item{dt}{(optional) the time step}
}
\value{
the wavelet transform
}
\description{
Perform Morlet Wavelet transform
}
