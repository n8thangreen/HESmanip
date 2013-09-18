\name{plot.FollowUpChart}
\alias{plot.FollowUpChart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.FollowUpChart}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.FollowUpChart(id, start, stop, col = 1, lty = 1, lwd = 2, xlim, xlab = "", ylab = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{stop}{
%%     ~~Describe \code{stop} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
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
\code{\link{plot.FollowUpChartBatch}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (id, start, stop, col = 1, lty = 1, lwd = 2, xlim, xlab = "", 
    ylab = "", ...) 
{
    if (missing(xlim)) 
        xlim <- c(min(start), max(stop))
    N <- length(start)
    plot(0, 0, type = "n", ylim = c(1, N/2), xlim = xlim, xlab = "", 
        ylab = "", ...)
    nix <- lapply(1:N, function(i) {
        segments(0, id[i], start[i], id[i], col = "lightgray", 
            lty = lty, lwd = lwd)
        segments(start[i], id[i], stop[i], id[i], col = col, 
            lty = lty, lwd = lwd)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
