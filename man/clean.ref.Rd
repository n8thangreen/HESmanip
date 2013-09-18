\name{clean.ref}
\alias{clean.ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clean reference table codes}
\description{
Clean reference table codes, removing extra characters and spaces.
}
\usage{
clean.ref(refcodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refcodes}{
%%     ~~Describe \code{refcodes} here~~
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
function (refcodes) 
{
    refcodes <- as.character(refcodes)
    refcodes <- gsub("'", "", refcodes)
    refcodes <- gsub(" ", "", refcodes)
    refcodes <- gsub(" ", "", refcodes)
    refcodes <- gsub(" ", "", refcodes)
    refcodes <- gsub("+", "", refcodes, fixed = TRUE)
    refcodes <- gsub("*", "", refcodes, fixed = TRUE)
    refcodes <- gsub(".", "", refcodes, fixed = TRUE)
    refcodes <- gsub("-", "", refcodes, fixed = TRUE)
    refcodes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
