\name{table.HR}
\alias{table.HR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{table.HR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
table.HR(output.HR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.HR}{
%%     ~~Describe \code{output.HR} here~~
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
\code{\link{table2.HR}},\code{\link{table3.HR}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (output.HR) 
{
    namesOrganisms <- names(output.HR)
    namesMethods <- names(output.HR[[1]])
    colNames <- c("organism", "method", "type", "beta", "exp(beta)", 
        "Lower CI", "Upper CI", "p")
    table.HR <- data.frame(matrix(ncol = length(colNames)))
    for (org in namesOrganisms) {
        for (meth in namesMethods) {
            namesEvent <- names(output.HR[[org]][[meth]])
            for (event in namesEvent) {
                table.HR <- rbind(table.HR, c(org, meth, event, 
                  extractCox(output.HR[[org]][[meth]][[event]])))
            }
        }
    }
    colnames(table.HR) <- colNames
    table.HR <- table.HR[!is.na(table.HR[, 1]), ]
    table.HR
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
