\name{HCAIsummary}
\alias{HCAIsummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HCAIsummary}
\description{
Aggregate the operations for summary statistics for age, length of stay & BSI time. Generally unformatted for reports.
}
\usage{
HCAIsummary(survData)
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
function (survData) 
{
    trimPerc <- 0.1
    print(c("age", summary(survData$age)))
    print(sd(survData$age, na.rm = TRUE))
    print(table(survData$age))
    print(c("LOS", summary(survData$time)))
    print(c("sd", sd(survData$time, na.rm = TRUE)))
    print(c("trimmed mean", mean(survData$time, trim = trimPerc)))
    print(c("spectime", summary(survData$spectime)))
    print(c("sd", sd(survData$spectime, na.rm = TRUE)))
    print(c("deaths", table(survData$event)))
    print(c("gender", table(survData$gender)))
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
