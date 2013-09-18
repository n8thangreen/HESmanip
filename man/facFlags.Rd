\name{facFlags}
\alias{facFlags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Group patients across multiple columns}
\description{
Group patients against specified groupings using change in factor labels, over all groups.
}
\usage{
facFlags(lookup, scheme, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lookup}{
Look-up table of codes by group.
}
  \item{scheme}{
Denotes which type of codes to use e.g. ICD10 or OPCS.
}
  \item{data}{
Patient data on which to apply the grouping.
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
\code{\link{my_recode}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lookup, scheme, data) 
{
    lookup.list <- as.list(unstack(lookup[lookup$Scheme == scheme, 
        c("Code", "Group")], Code ~ Group))
    y <- apply(data, 2, function(x) my_recode(x, lookup.list))
    list.names <- names(lookup.list)
    res <- NULL
    for (i in list.names) {
        res <- cbind(res, apply(y, 1, function(x) any(x == i)))
    }
    colnames(res) <- list.names
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
