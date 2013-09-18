\name{commonlist}
\alias{commonlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{commonlist}
\description{
Produces an ordered list by frequency.
}
\usage{
commonlist(group, cut = 20, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{cut}{
frequency cut-off value
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (group, cut = 20, plot = FALSE) 
{
    x <- as.data.frame(table(group))
    x <- data.frame(x, prop = x$Freq/sum(x$Freq))
    if (plot) {
        plot(x$Var1[x$Freq > cut], x$Freq[x$Freq > cut], las = 3)
    }
    x[x$Freq > cut, ][order(x[x$Freq > cut, "Freq"], decreasing = TRUE), 
        ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
