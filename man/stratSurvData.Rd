\name{stratSurvData}
\alias{stratSurvData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{stratSurvData}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stratSurvData(survData, survDataCol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{survDataCol}{
%%     ~~Describe \code{survDataCol} here~~
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
function (survData, survDataCol) 
{
    Glevels <- nlevels(survDataCol)
    Gnames <- levels(survDataCol)
    survDataStrat <- list()
    for (i in 1:Glevels) {
        survDataStrat[[i]] <- survData[survDataCol == Gnames[i], 
            ]
        survDataStrat[[i]] <- survDataStrat[[i]][order(survDataStrat[[i]]$infstatus), 
            ]
    }
    names(survDataStrat) <- Gnames
    survDataStrat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
