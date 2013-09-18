\name{boot.clos.dtable}
\alias{boot.clos.dtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap LOS with data.table}
\description{
A data.table version of boot_clos() for faster run-time.
}
\usage{
boot.clos.dtable(data, state.names, tra, cens.name, s = 0, nboot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  class(data.table)
}
  \item{state.names}{
%%     ~~Describe \code{state.names} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{cens.name}{
%%     ~~Describe \code{cens.name} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
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
\code{\link{boot.clos}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, state.names, tra, cens.name, s = 0, nboot) 
{
    res <- double(nboot)
    setkey(data, "id")
    for (i in seq_len(nboot)) {
        index <- sample(unique(data[, id]), replace = TRUE)
        dboot <- data[J(index)]
        tr.prob <- etm(dboot, state.names, tra, cens.name, s, 
            cova = FALSE)
        res[i] <- etm::clos(tr.prob)$e.phi
    }
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
