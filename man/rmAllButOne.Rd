\name{rmAllButOne}
\alias{rmAllButOne}
\title{Remove all but one patient record}
\description{
Remove all but one patient record to omit duplicates.
}
\usage{
rmAllButOne(ind, subrows, total.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
%%     ~~Describe \code{ind} here~~
}
  \item{subrows}{
row names of subset of rows of \code{total.data} on which to apply function.
}
  \item{total.data}{
Master data.frame of patient records
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ind, subrows, total.data) 
{
    ind.min <- min(which(ind[subrows]))
    ind[subrows] <- FALSE
    ind[subrows[ind.min]] <- TRUE
    total.data <- total.data[ind, ]
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
