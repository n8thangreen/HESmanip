\name{listLOS}
\alias{listLOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of stay multistate model statistics}
\description{
Calculate length of stay multistate model statistics.
}
\usage{
listLOS(survData, tra, type, standerr = FALSE, trim = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{tra}{
Transition matrix of permitted direct between state jumps.
}
  \item{type}{
Model type, either cause-specific, subsdistribution with discharge alive as event of interest or in-hospital death as event of interest.
}
  \item{standerr}{
Logical, whether the standard error should be calculated.
}
  \item{trim}{
Upper percentile of data to trim.
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
function (survData, tra, type, standerr = FALSE, trim = 0.99) 
{
    tsurvData <- survData
    msm.data <- etmFormat(tsurvData, type)
    mvna.data <- mvna::mvna(data = msm.data, state.names = c("0", 
        "1", "2"), tra = tra, cens.name = "cens")
    etm.data <- etm::etm(data = msm.data, state.names = c("0", 
        "1", "2"), tra = tra, cens.name = "cens", s = 0)
    clos.data <- etm::clos(etm.data, aw = TRUE)
    if (standerr) {
        se <- sqrt(var(boot.clos(data = msm.data, state.names = c("0", 
            "1", "2"), tra = tra, cens.name = "cens", 0, nboot = 10)))
    }
    else {
        se <- 0
    }
    list(mvna.data = mvna.data, clos.data = clos.data, se = se, 
        etm.data = etm.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
