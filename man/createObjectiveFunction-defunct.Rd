% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrar-DEPRECATED.R
\name{createObjectiveFunction-defunct}
\alias{createObjectiveFunction-defunct}
\title{Create an objective function to be used with optimization routines}
\arguments{
\item{runModel}{Function to run the model and produce a list of outputs.}

\item{info}{A data.frame with the information about the calibration, 
normally created with the \code{\link{getCalibrationInfo}} function. 
See details.}

\item{observed}{A list of the observed variables created with the 
function \code{\link{getObservedData}}}

\item{aggFn}{A function to aggregate \code{fn} to a scalar value if the
returned value is a vector. Some optimization algorithm can explote the
additional information provided by a vectorial output from \code{fn}}

\item{aggregate}{boolean, if TRUE, a scalar value is returned using the 
\code{aggFn}.}

\item{\dots}{More arguments passed to the \code{runModel} function.}
}
\value{
A function, integrating the simulation of the model and the 
comparison with observed data.
}
\description{
Create a new function, to be used as the objective function in the 
calibration, given a function to run the model within R, observed data 
and information about the comparison with data.
}
\seealso{
\code{\link{calibrar-defunct}}
}
\author{
Ricardo Oliveros-Ramos
}
\keyword{internal}
