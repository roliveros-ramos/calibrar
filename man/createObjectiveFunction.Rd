\name{createObjectiveFunction}
\alias{createObjectiveFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createObjectiveFunction(runModel, info, observed, aggFn = .weighted.sum, aggregate = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{runModel}{
%%     ~~Describe \code{runModel} here~~
}
  \item{info}{
%%     ~~Describe \code{info} here~~
}
  \item{observed}{
%%     ~~Describe \code{observed} here~~
}
  \item{aggFn}{
%%     ~~Describe \code{aggFn} here~~
}
  \item{aggregate}{
%%     ~~Describe \code{aggregate} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
function (runModel, info, observed, aggFn = .weighted.sum, aggregate = FALSE, 
    ...) 
{
    fn = match.fun(runModel)
    aggFn = match.fun(aggFn)
    force(observed)
    force(info)
    force(aggregate)
    fn1 = function(par) {
        aggFn = match.fun(aggFn)
        simulated = fn(par, ...)
        output = .calculateObjetiveValue(obs = observed, sim = simulated, 
            info = info)
        if (isTRUE(aggregate)) 
            output = aggFn(x = output, w = info$weights)
        return(output)
    }
    return(fn1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
