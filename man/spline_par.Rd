% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrar-splines.R
\name{spline_par}
\alias{spline_par}
\title{Predict time-varying parameters using splines.}
\usage{
spline_par(par, n, knots = NULL, periodic = FALSE, period = NULL)
}
\arguments{
\item{par}{Values at knots}

\item{n}{Number of points. Time (independent variable) is assumed to be between 0 and n with length(par) equidistant points (including 0 and n).}

\item{knots}{Position of knots. Default, is length(x) equidistant points between 0 and 1. Always are re-scaled to 0 to 1.}

\item{periodic}{boolean, is the spline periodic?}

\item{period}{If periodic is TRUE, it specify the time period.}
}
\value{
A list with the interpolates values as 'x' and 'time'.
}
\description{
Predict time-varying parameters using splines.
}
