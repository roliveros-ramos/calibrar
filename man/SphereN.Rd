\name{SphereN}
\alias{SphereN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sphere function with random noise
%%  ~~function to do ... ~~
}
\description{
This function calculates the Euclidian distance from a point to the origin after a random displacement of it position.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SphereN(x, sd = 0.1, aggregate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
The coordinates of the point
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
The standard deviation of the noise to be added to the position of \code{x},
a normal distribution with mean zero is used.
}
  \item{aggregate}{
%%     ~~Describe \code{aggregate} here~~
If \code{aggregate} is \code{TRUE} the distance is returned, otherwise the size of the projection of the distance among each axis.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
The distance from the point \code{x} to the origin after a random displacement.
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Ricardo Oliveros--Ramos
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, sd = 0.1, aggregate = TRUE) 
{
    x = x + rnorm(length(x), sd = sd)
    out = x^2
    if (isTRUE(aggregate)) 
        return(sum(out))
    else return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{stochastic}
\keyword{random}% __ONLY ONE__ keyword per line
