\name{hr.tv.causespecific}
\alias{hr.tv.causespecific}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time-dependent Cox regression for the cause-specific model}
\description{
Time-dependent Cox regression for the cause-specific model.
}
\usage{
hr.tv.causespecific(survData)
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
\code{\link{hr.naive}}, \code{\link{hr.tv.subdistribution}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (survData) 
{
    tvdata <- tvCox(survData, type = "")
    y.alive <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status == 
        1)
    fit.alive.timeonly <- summary(coxph(y.alive ~ cluster(id) + 
        inf, data = tvdata))
    fit.alive.full <- summary(coxph(y.alive ~ cluster(id) + age + 
        sex + Tai + other + head + cath + surgical + strata(prem) + 
        cancer + hes_neocaredescription + hes_admimethdescription + 
        hes_admisorcdescription + inf, data = tvdata))
    y.dead <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status == 
        2)
    fit.dead.timeonly <- summary(coxph(y.dead ~ cluster(id) + 
        inf, data = tvdata))
    fit.dead.full <- summary(coxph(y.dead ~ cluster(id) + age + 
        sex + Tai + other + head + cath + surgical + strata(prem) + 
        cancer + hes_neocaredescription + hes_admimethdescription + 
        hes_admisorcdescription + inf, data = tvdata))
    y.both <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status != 
        0)
    fit.both.timeonly <- summary(coxph(y.both ~ cluster(id) + 
        inf, data = tvdata))
    fit.both.full <- summary(coxph(y.both ~ cluster(id) + age + 
        sex + Tai + other + head + cath + surgical + strata(prem) + 
        cancer + hes_neocaredescription + hes_admimethdescription + 
        hes_admisorcdescription + inf, data = tvdata))
    return(list(fit.alive.timeonly = fit.alive.timeonly, fit.alive.full = fit.alive.full, 
        fit.dead.timeonly = fit.dead.timeonly, fit.dead.full = fit.dead.full, 
        fit.both.timeonly = fit.both.timeonly, fit.both.full = fit.both.full))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
