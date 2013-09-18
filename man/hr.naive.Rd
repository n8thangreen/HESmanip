\name{hr.naive}
\alias{hr.naive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time-independent Cox regression}
\description{
Time-independent Cox regression for alternative models.
}
\usage{
hr.naive(survData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{fit.disch}{Discharge alive as event of interest.}
  \item{fit.dead}{In-hospital death as event of interest.}
  \item{fit.both}{Combined end-point of discharged alive and in-hospital death.}
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
\code{\link{hr.tv.causespecific}}, \code{\link{hr.tv.subdistribution}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (survData) 
{
    survData$spectime <- ifelse(is.na(survData$spectime), 0, 
        survData$spectime)
    cov <- paste(c(names(survData)[1:17], "age", "gender", "infstatus"), 
        collapse = "+")
    fit.disch <- summary(coxph(Surv(time, 1 - as.logical(event)) ~ 
        age + gender + Tai + other + head + cath + surgical + 
            prem + cancer + hes_neocaredescription + hes_admimethdescription + 
            hes_admisorcdescription + infstatus, data = survData))
    survData$event[survData$missingType == 0] <- FALSE
    fit.dead <- summary(coxph(Surv(time, as.logical(event)) ~ 
        age + gender + Tai + other + head + cath + surgical + 
            prem + cancer + hes_neocaredescription + hes_admimethdescription + 
            hes_admisorcdescription + infstatus, data = survData))
    fit.both <- summary(coxph(Surv(time, missingType) ~ age + 
        gender + Tai + other + head + cath + surgical + prem + 
        cancer + hes_neocaredescription + hes_admimethdescription + 
        hes_admisorcdescription + infstatus, data = survData))
    return(list(fit.disch = fit.disch, fit.dead = fit.dead, fit.both = fit.both))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
