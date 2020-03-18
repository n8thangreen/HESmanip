\name{ICDwordmatch}
\alias{ICDwordmatch}
\title{ICD10 word match against patient codes}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ICDwordmatch(ICD10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ICD10}{
Array with \code{Group}, \code{Code} and \code{Description} columns.
}
}
\details{
Regular expression to search in Description and Group names for a word match to the hard-coded risk factor list, and then return the relevant Code. The risk factor words to match are cancer, cathet, dialys, low birth, preterm, prematurity, intravenous.
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
function (ICD10) 
{
    RFcodes <- list()
    cat <- c("cancer", "cathet", "dialys", "low birth", "preterm", 
        "prematurity", "^intravenous")
    for (i in seq_along(cat)) {
        catmod <- paste("(?i)", cat[i], sep = "")
        RFrows <- unique(c(grep(catmod, ICD10$Group, fixed = FALSE), 
            grep(catmod, ICD10$Description, fixed = FALSE)))
        RFcodes[[i]] <- ICD10$Code[RFrows]
    }
    names(RFcodes) <- cat
    RFcodes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
