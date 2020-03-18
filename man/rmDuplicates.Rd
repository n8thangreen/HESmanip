\name{rmDuplicates}
\alias{rmDuplicates}
\title{Remove duplicate patient records}
\description{
Remove patient records for the same patient-stay where one or more of the records have duplicate or missing entries.
}
\usage{
rmDuplicates(total.data, subrows)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
Complete combined patient record data set.
}
  \item{subrows}{
Row names of subset of patient records on which to apply function.
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
\code{\link{findDuplicates}}
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
