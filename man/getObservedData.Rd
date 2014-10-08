\name{getObservedData}
\alias{getObservedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getObservedData(info, path, data.folder = "data", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{info}{
%%     ~~Describe \code{info} here~~
}
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{data.folder}{
%%     ~~Describe \code{data.folder} here~~
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
function (info, path, data.folder = "data", ...) 
{
    observed = list()
    variables = info$variable
    useData = as.logical(info$useData)
    cat("Creating observed data list for calibration...", "\n")
    for (var in 1:nrow(info)) {
        cat(paste0("Variable: ", variables[var], "\n"))
        var.path = file.path(path, data.folder, paste0(variables[var], 
            ".csv"))
        datos = if (useData[var]) 
            .read.csv3(var.path, ...)
        else NA
        observed[[var]] = datos
    }
    names(observed) = variables
    return(observed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
