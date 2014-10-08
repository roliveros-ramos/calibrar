\name{getCalibrationInfo}
\alias{getCalibrationInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCalibrationInfo(path, file = "calibrationInfo.csv", stringsAsFactors = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{stringsAsFactors}{
%%     ~~Describe \code{stringsAsFactors} here~~
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
function (path, file = "calibrationInfo.csv", stringsAsFactors = FALSE, 
    ...) 
{
    caliPath = file.path(path, file)
    calibrationInfo = read.csv(caliPath, stringsAsFactors = FALSE, 
        ...)
    fullNames = c("variable", "type", "calibrate", "weights", 
        "useData")
    doesNotMatch = !(names(calibrationInfo) \%in\% fullNames)
    dnm = names(calibrationInfo)[doesNotMatch]
    isMissing = !(fullNames \%in\% names(calibrationInfo))
    im = fullNames[isMissing]
    sdnm = if (length(dnm) > 1) 
        " columns do "
    else " column does "
    sim = if (length(im) > 1) 
        " variables are "
    else " variable is "
    msg1 = paste0("Error in ", caliPath, " file (", paste(sapply(dnm, 
        sQuote), collapse = ", "), sdnm, "not match).")
    msg2 = paste0("Error in ", caliPath, " file (", paste(sapply(im, 
        sQuote), collapse = ", "), sim, "missing).")
    if (any(doesNotMatch)) 
        stop(msg1)
    if (any(isMissing)) 
        stop(msg2)
    calibrationInfo$variable = as.character(calibrationInfo$variable)
    calibrationInfo$type = as.character(calibrationInfo$type)
    calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
    calibrationInfo$weights = as.numeric(calibrationInfo$weights)
    calibrationInfo$useData = as.logical(calibrationInfo$useData)
    return(calibrationInfo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
