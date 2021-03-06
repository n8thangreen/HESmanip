\name{subdistn}
\alias{subdistn}
\title{Fine & Gray subdistribution model time to event times}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subdistn(survData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
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
A Proportional Hazards Model for the Subdistribution of a Competing Risk, Jason P. Fine and Robert J. Gray, Journal of the American Statistical Association, Vol. 94, No. 446 (Jun 1999), pp. 496-509.
}
\author{
N Green
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
function (survData) 
{
    survData$time2disch <- survData$time2death <- survData$time
    mtime1 <- max(survData[survData$infstatus == 1, "time"], 
        na.rm = TRUE)
    survData$time2disch[survData$infstatus == 1 & survData$event & 
        survData$missingType] <- mtime1 + 1
    survData$time2death[survData$infstatus == 1 & !survData$event & 
        survData$missingType] <- mtime1 + 1
        
    mtime0 <- max(survData[survData$infstatus == 0, "time"], 
        na.rm = TRUE)
    survData$time2disch[survData$infstatus == 0 & survData$event & 
        survData$missingType] <- mtime0 + 1
    survData$time2death[survData$infstatus == 0 & !survData$event & 
        survData$missingType] <- mtime0 + 1
    survData
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
