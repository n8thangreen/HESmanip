\name{highlevelHCAI}
\alias{highlevelHCAI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combined length of stay calculations}
\description{
Combine length of stay calculation and Cox proportional hazard ratio calculation across defined groups.
}
\usage{
highlevelHCAI(survData.list)
}
\arguments{
  \item{survData.list}{A list of patient record data stratified by some comorbidity or confounder of interest.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{output.HR}{Cox regression output for the all model types.}
  \item{ouput.LOS}{Length of stay output using multistate model.}
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (survData.list) 
{
    output.HR <- byOrganismHR(survData.list)
    type <- ""
    survData.list <- lapply(survData.list, function(x) x[!is.na(x$time) & 
        x$time > 0, ])
    survData.list.rm <- lapply(survData.list, function(x) x[x$time >= 
        min(x$spectime[x$infstatus == 1]), ])
    output.LOS <- bygroupLOS(survData.list.rm, type, standerr = TRUE)
    list(output.HR = output.HR, output.LOS = output.LOS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
