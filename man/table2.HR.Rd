\name{table2.HR}
\alias{table2.HR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{table2.HR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
table2.HR(res, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
\code{\link{table.HR}},\code{\link{table3.HR}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res, model) 
{
    namesGroup <- unique(res$organism)
    numGroup <- length(namesGroup)
    res.sub <- res[res$method == model, c("organism", "exp(beta)", 
        "Lower CI", "Upper CI")]
    subHeads <- c("HR", "LCI", "UCI")
    colHeads <- c("organism", paste("atime", subHeads), paste("afull", 
        subHeads), paste("dtime", subHeads), paste("dfull", subHeads))
    res.new <- data.frame(matrix(ncol = length(colHeads), nrow = numGroup), 
        check.rows = FALSE)
    names(res.new) <- colHeads
    for (j in 1:numGroup) {
        res.temp <- NA
        firstrow <- min(which(res.sub$organism == namesGroup[j]))
        for (i in 1:4) {
            res.temp <- cbind(res.temp, res.sub[firstrow + i - 
                1, -1])
        }
        res.new[j, ] <- res.temp
    }
    res.new[, 1] <- namesGroup
    res.new
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
