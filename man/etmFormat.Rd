\name{etmFormat}
\alias{etmFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform to etm format}
\description{
Transform data.frame in to package:etm \code{etm} (long) format.
}
\usage{
etmFormat(survData, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{type}{
Model type.
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
function (survData, type) 
{
    tvdata <- tvCox(survData, type = type)
    tvdata <- tvdata[order(tvdata$originf, decreasing = TRUE), 
        ]
    numCases.inf <- sum(tvdata$originf == 1)
    numCases.mix <- sum(tvdata$originf == 0)
    id <- tvdata$id
    entry <- tvdata$tstart
    exit <- tvdata$tstop
    from <- c(rep(0:1, length.out = numCases.inf), rep(0, numCases.mix))
    to <- c(rep(1:2, length.out = numCases.inf), rep(2, numCases.mix))
    to2 <- tvdata$disch + to
    rinf <- 1:numCases.inf
    rmix <- (max(rinf) + 1):nrow(tvdata)
    rcens.inf <- tvdata$inf[rinf] == 1 & (tvdata$disch[rinf] + 
        tvdata$death[rinf]) == 0
    rcens.mix <- (tvdata$disch[rmix] + tvdata$death[rmix]) == 
        0
    rcens <- c(rcens.inf, rcens.mix)
    to[rcens] <- to2[rcens] <- "cens"
    age <- tvdata$age
    agegr <- tvdata$agegr
    sex <- tvdata$sex
    msm.data <- data.frame(id, from, to, entry, exit)
    msm.data <- msm.data[!is.na(msm.data$exit), ]
    msm.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
