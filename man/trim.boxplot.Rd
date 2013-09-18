\name{trim.boxplot}
\alias{trim.boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{trim.boxplot}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trim.boxplot(survData, trim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{trim}{
%%     ~~Describe \code{trim} here~~
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
function (survData, trim) 
{
    for (i in 1:length(trim)) {
        out.trim[[i]] <- listLOS(survData, tra, type, standerr = TRUE, 
            trim = trim[i])
        names(out.trim) <- trim
    }
    byGroupLOS.boxplot(out.trim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
