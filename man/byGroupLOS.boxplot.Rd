\name{byGroupLOS.boxplot}
\alias{byGroupLOS.boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot of length of stay by groups}
\description{
Box plots for length of stay for each specified group with error bars.
}
\usage{
byGroupLOS.boxplot(output.LOS, offset = 0, ADD = FALSE, axis = "s", col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.LOS}{
%%     ~~Describe \code{output.LOS} here~~
}
  \item{offset}{
How much to offset the bars from integer values. Used when plotting multiple bars on a single plot.
}
  \item{ADD}{
Indicates if this plot is to be overlayed on to a previous plot.
}
  \item{axis}{
Default is to label axis but this can be suppressed with "n".
}
  \item{col}{
Colour of bars.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NULL
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
\code{\link{byGroupLOS}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (output.LOS, offset = 0, ADD = FALSE, axis = "s", col = "black") 
{
    namesGroup.short <- names(output.LOS)
    wmeans <- se <- NA
    numGroup <- length(output.LOS)
    neworder <- 1:numGroup
    for (i in 1:numGroup) {
        if (is.null(output.LOS[[i]]$clos.data$e.phi)) {
            wmeans[i] <- se[i] <- 0
        }
        else {
            wmeans[i] <- output.LOS[[i]]$clos.data$e.phi
            se[i] <- output.LOS[[i]]$se
        }
    }
    wmeans <- unlist(wmeans)[neworder]
    se <- se[neworder]
    par(mar = c(5, 9, 2, 0.5))
    bp <- boxplot(0, plot = FALSE, color = "yellow")
    bp$stats <- rbind(wmeans - (2.33 * se), wmeans - (1.96 * 
        se), wmeans, wmeans + (1.96 * se), wmeans + (2.33 * se))
    rownames(bp$stats) <- c("", "", "", "", "")
    bp$n <- rep(20, numGroup)
    bp$conf <- rbind(wmeans - (2.33 * se), wmeans + (2.33 * se))
    bp$names <- namesGroup.short[neworder]
    bxp(bp, horizontal = T, xlab = "Excess length of stay (days)", 
        boxwex = 0.25, las = 1, add = ADD, at = (1:numGroup) + 
            offset, border = col, yaxt = axis)
    abline(v = 0, col = "red", lty = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
