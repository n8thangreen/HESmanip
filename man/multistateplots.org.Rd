\name{multistateplots.org}
\alias{multistateplots.org}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{multistateplots.org}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multistateplots.org(output.LOS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.LOS}{
%%     ~~Describe \code{output.LOS} here~~
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
function (output.LOS) 
{
    for (i in names(output.LOS)) {
        plot(output.LOS[[i]]$clos.data, xlim = c(0, 40), ylim.e = c(0, 
            90), main = i)
        plot(output.LOS[[i]]$mvna.data, xlim = c(0, 90), main = i)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
