\name{plot.FollowUpChartBatch}
\alias{plot.FollowUpChartBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.FollowUpChartBatch}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.FollowUpChartBatch(survData.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData.list}{
%%     ~~Describe \code{survData.list} here~~
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
\code{\link{plot.FollowUpChart}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (survData.list) 
{
    namesGroup <- c("All", "Gram-positive", "Gram-negative", 
        "CoNS", expression(paste(italic("Enterococcus"), plain(spp.))), 
        expression(italic("S. aureus")), "Other", expression(italic("E. Coli")), 
        "Non-p. Streptococci", expression(paste(italic("Klebsiella"), 
            plain(spp.))), expression(paste(italic("Enterobacter"), 
            plain(spp.))))
    numgroups <- c(1, 4, 5, 6, 8, 9, 10, 11)
    par(mar = c(5, 6, 2, 1))
    par(oma = c(3, 3, 3, 0))
    par(mai = c(0.5, 0.6, 0.3, 0.3))
    par(mfrow = c(4, 2))
    for (i in numgroups) {
        data <- tvCox(survData.list[[i]][order(survData.list[[i]]$time), 
            ], type = "")
        data <- data[data$originf == 1, ]
        plot.FollowUpChart(data$id, start = data$tstart, stop = data$tstop, 
            xlab = "Time since admission (days)", ylab = "HA-BSI patients", 
            main = namesGroup[i], xlim = c(0, 250))
        points(data$tstop, data$id, pch = "l")
        points(data$tstop[data$disch == 1], data$id[data$disch == 
            1], pch = 19, cex = 1)
        points(data$tstop[data$death == 1], data$id[data$death == 
            1], pch = 8, cex = 1, lwd = 2, col = "red")
    }
    mtext("Time since admission (days)", side = 1, outer = T)
    mtext("HA-BSI patients", side = 2, outer = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
