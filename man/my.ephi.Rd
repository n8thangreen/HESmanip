\name{my.ephi}
\alias{my.ephi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{my.ephi}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my.ephi(clos.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clos.data}{
%%     ~~Describe \code{clos.data} here~~
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
function (clos.data) 
{
    clos.new <- merge(cbind(clos.data$time, clos.data$phi.case, 
        clos.data$phi.control), cbind(clos.data$w.time, clos.data$weights), 
        by = 1)
    return((clos.new[, 2] - clos.new[, 3]) \%*\% clos.new[, 4])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
