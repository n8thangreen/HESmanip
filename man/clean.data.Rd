\name{clean.data}
\alias{clean.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clean ICD10 codes}
\description{
In ICD10 codes remove unwanted trailing characters.
}
\usage{
clean.data(data, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
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
function (data, cols) 
{
    for (i in cols) {
        data[, i] <- gsub("-$", "0", data[, i])
        data[, i] <- gsub("X$", "0", data[, i])
        data[, i] <- gsub("X.$", "0", data[, i])
        data[, i] <- substring(data[, i], 1, 4)
    }
    data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
