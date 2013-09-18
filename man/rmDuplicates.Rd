\name{rmDuplicates}
\alias{rmDuplicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rmDuplicates}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmDuplicates(total.data, subrows)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
%%     ~~Describe \code{total.data} here~~
}
  \item{subrows}{
%%     ~~Describe \code{subrows} here~~
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
function (total.data, subrows) 
{
    inc.ind <- hos.ind <- rep(TRUE, nrow(total.data))
    inc.ind[subrows] <- total.data$hes_dischargeddeadoralive[subrows] \%in\% 
        c("Death", "Live")
    hos.ind[subrows] <- total.data$hes_dismethdescription[subrows] \%in\% 
        c("Not applicable patient still in hospital", "Not known")
    if (any(inc.ind[subrows])) {
        total.data <- rmAllButOne(inc.ind, subrows, total.data)
    }
    else if (all(hos.ind[subrows])) {
        if (any(!is.na(total.data$hes_disdte[subrows]))) {
            hos.ind[subrows] <- total.data$hes_disdt[subrows] == 
                max(total.data$hes_disdte[subrows], na.rm = TRUE)
        }
        total.data <- rmAllButOne(hos.ind, subrows, total.data)
    }
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
