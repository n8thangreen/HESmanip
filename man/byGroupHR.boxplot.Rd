\name{byGroupHR.boxplot}
\alias{byGroupHR.boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hazard ratio boxplots}
\description{
Boxplot of hazard ratio by specified groups.
}
\usage{
byGroupHR.boxplot(res2, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res2}{
  Of format table2()
}
  \item{model}{
Particular model type; either naive, subdistribution or cause-specific.
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
N Green
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{HRboxplot.batch}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res2, model) 
{
    namesGroup.short <- c("CoNS", "Other", "Non-p. Streptococci", 
        expression(italic("S. aureus")), expression(paste(italic("Enterococcus"), 
            plain(spp.))), expression(italic("E. Coli")), expression(paste(italic("Klebsiella"), 
            plain(spp.))), expression(paste(italic("Enterobacter"), 
            plain(spp.))), "Gram-negative", "Gram-positive", 
        "All")
    col <- grep(model, names(res2))
    if (model == "atime HR") {
        model <- "Discharge with time to infection only"
    }
    if (model == "afull HR") {
        model <- "Discharge with full set of covariates"
    }
    if (model == "dtime HR") {
        model <- "Death with time to infection only"
    }
    if (model == "dfull HR") {
        model <- "Death with full set of covariates"
    }
    model <- ""
    neworder <- c(11, 10, 9, 1, 6, 8, 5, 7, 3, 4, 2)
    wmeans <- as.numeric(res2[, col])[neworder]
    liw <- as.numeric(res2[, col + 1])[neworder]
    uiw <- as.numeric(res2[, col + 2])[neworder]
    par(mar = c(10, 4, 2, 0.5))
    plotCI(barplot(wmeans, col = "blue", ylim = c(0, min(10, 
        max(uiw, rm.na = TRUE))), names = namesGroup.short[neworder], 
        las = 2, main = model, ylab = "Hazard ratio"), wmeans, 
        uiw = uiw - wmeans, liw = wmeans - liw, add = TRUE)
    abline(h = 1, lty = 2, col = "red")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
