% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wavelet_morlet.R
\name{mo}
\alias{mo}
\title{Definition of the Morlet wavelet}
\usage{
mo(t, t_0 = 0, omega = 6, j = 0)
}
\arguments{
\item{t}{a numeric vector corresponding to times
where to evaluate the Morlet wavelet}

\item{t_0}{(optional) a time indicating the timing translation}

\item{j}{(optional) an index corresponding to the calculated scale s_j
(the corresponding scale is s_0  2^(j*dj))}
}
\value{
evaluation of the Morlet wavelet at values in \code{t}
}
\description{
Definition of the Morlet wavelet
}
\examples{
N_4_original(2)
}
