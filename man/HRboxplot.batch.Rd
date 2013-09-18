\name{HRboxplot.batch}
\alias{HRboxplot.batch}
\title{Hazard ratio boxplots for all groups as a batch run}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HRboxplot.batch(res2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res2}{
Array of type output \code{table2}.
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
\code{\link{byGroupHR.boxplot}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res2) 
{
    HRcol <- grep("HR", names(res2))
    HRnames <- names(res2)[HRcol]
    for (i in HRnames) {
        byGroupHR.boxplot(res2, i)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
