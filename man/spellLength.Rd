\name{spellLength}
\alias{spellLength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{spellLength}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spellLength(adm, disch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adm}{
%%     ~~Describe \code{adm} here~~
}
  \item{disch}{
%%     ~~Describe \code{disch} here~~
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
function (adm, disch) 
{
    for (i in 1:length(adm)) {
        dupl <- adm[i] == disch
        disch[dupl] <- disch[i]
        if (any(dupl, na.rm = TRUE)) {
            adm[i] <- adm[dupl] <- min(adm[dupl], na.rm = TRUE)
        }
    }
    list(adm = adm, disch = disch)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
