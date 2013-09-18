\name{my_recode}
\alias{my_recode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Group together codes}
\description{
According to a look-up table, patients are grouped together.
}
\usage{
my_recode(codes, lookuplist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{codes}{
%%     ~~Describe \code{codes} here~~
}
  \item{lookuplist}{
%%     ~~Describe \code{lookuplist} here~~
}
}
\details{
Uses a neat trick of relabelling factors which is much simpler and quicker than stepping through each patients records.
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
\code{\link{facFlags}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (codes, lookuplist) 
{
    nfac <- factor(codes)
    lfac <- levels(nfac)
    othrlevs <- lfac[!lfac \%in\% unlist(lookuplist)]
    levels(nfac) <- c(lookuplist, list(all_others = othrlevs))
    nfac
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
