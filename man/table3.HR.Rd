\name{table3.HR}
\alias{table3.HR}
\title{More condensed table format of hazard ratio values for different discharge states}
\description{
One of several output formats. This format is used in the Journal of Pediatric Infectious Diseases Society (JPIDS) paper.
}
\usage{
table3.HR(res, hrtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{hrtype}{
%%     ~~Describe \code{hrtype} here~~
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
\code{\link{table.HR}},\code{\link{table2.HR}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res, hrtype) 
{
    if (hrtype == "timedeptsubdistn") {
        colNames <- c("Organism", "Disch Time-adjusted", "Disch Fully adjusted", 
            "Death Time-adjusted", "Death Fully adjusted")
    }
    else if (hrtype == "naive") {
        colNames <- c("Organism", "Disch", "Death", "Both")
    }
    else if (hrtype == "timedeptcausespec") {
        colNames <- c("Organism", "Disch Time-adjusted", "Disch Fully adjusted", 
            "Death Time-adjusted", "Death Fully adjusted", "Both Time-adjusted", 
            "Both Fully adjusted")
    }
    else {
        stop("type unidentified")
    }
    res.new <- data.frame(matrix(ncol = length(colNames)))
    colnames(res.new) <- colNames
    for (name in unique(res$organism)) {
        whichrows <- which(res$organism == name & res$method == 
            hrtype)
        rowTotal <- NULL
        for (j in 1:length(whichrows)) {
            temp <- paste(res[whichrows[j], "exp(beta)"], " (", 
                res[whichrows[j], "Lower CI"], ", ", res[whichrows[j], 
                  "Upper CI"], ")", sep = "")
            rowTotal <- c(rowTotal, temp)
        }
        res.new <- rbind(res.new, c(name, rowTotal))
    }
    res.new <- res.new[!is.na(res.new[, 1]), ]
    res.new
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
