\name{boot.clos}
\alias{boot.clos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap LOS}
\description{Bootstrap standard error estimate for length of stay.}
\usage{
boot.clos(data, state.names, tra, cens.name, s = 0, nboot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
\code{\link{boot.clos.dtable}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, state.names, tra, cens.name, s = 0, nboot) 
{
    res <- double(nboot)
    for (i in seq_len(nboot)) {
        index.inf <- sample(x = unique(data$id[data$to == 1]), 
            replace = TRUE)
        index.mix <- sample(x = unique(data$id[data$from == 0 & 
            data$to == 2]), replace = TRUE)
        index <- c(index.inf, index.mix)
        linds <- sapply(index, function(x) which(x == data$id))
        indrep <- sapply(linds, function(y) length(y))
        inds <- unlist(linds)
        new.id <- rep(seq_along(index), indrep)
        dboot <- cbind(data[inds, ], new.id)
        dboot$id <- dboot$new.id
        dboot$new.id <- NULL
        tr.prob <- etm::etm(dboot, state.names, tra, cens.name, 
            s, cova = FALSE)
        res[i] <- etm::clos(tr.prob)$e.phi
    }
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
