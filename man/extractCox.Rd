\name{extractCox}
\alias{extractCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract subset of Cox regression statistics.}
\description{
Extract a subset of the Cox proportional hazard output values and format in report-ready table.
}
\usage{
extractCox(cox)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cox}{
%%     ~~Describe \code{cox} here~~
}
}
\details{
The infection status variable assumed to be the last covariate in list.
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
function (cox) 
{
    dpl <- 3
    lastRow <- nrow(coef(cox))
    beta <- round(coef(cox)[lastRow, "coef"], dpl)
    se <- round(coef(cox)[lastRow, "se(coef)"], dpl)
    p <- round(coef(cox)[lastRow, "Pr(>|z|)"], dpl)
    CI <- round(cox$conf.int[lastRow, c("lower .95", "upper .95")], 
        dpl)
    res <- cbind(beta, `exp(beta)` = round(exp(beta), 3), CI[1], 
        CI[2], p)
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
