\name{trimSurv}
\alias{trimSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{trimSurv}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trimSurv(survData, percentile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{percentile}{
%%     ~~Describe \code{percentile} here~~
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
function (survData, percentile) 
{
    qdata.mix <- quantile(survData[survData$infstatus == 0, ]$time, 
        c(1 - percentile, percentile))
    qdata.inf <- quantile(survData[survData$infstatus == 1, ]$time, 
        c(1 - percentile, percentile))
    survData <- survData[(survData$time <= qdata.mix[2] & survData$infstatus == 
        0) | (survData$infstatus == 1), ]
    survData
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
