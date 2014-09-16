\name{calibraR-package}
\alias{calibraR-package}
\alias{calibraR}
\docType{package}
\title{
Calibration of Ecological Models using Evolutionary Algorithms
}
\description{
Calibration of Ecological Models using Evolutionary Algorithms
}
\details{
\tabular{ll}{
Package: \tab calibraR\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-09-14\cr
License: \tab GPL-2\cr
}
calibrate()
optimES()
}
\author{
Ricardo Oliveros-Ramos
Maintainer: Ricardo Oliveros-Ramos <ricardo.oliveros@gmail.com>
}
\references{
calibrar: an R package for the calibration of ecological models
}
\keyword{ calibration }

\examples{
calibrate(par=rep(NA, 5), fn=SphereN)
calibrate(par=rep(NA, 5), fn=SphereN, replicates=3)
calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5)
calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
calibrate(par=rep(0.5, 5), fn=SphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
}

