\name{byGroupHR}
\alias{byGroupHR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Cox regression by groups}
\description{
Batch calculate all the Cox regressions for each group in turn.
Can select which type of hazard ratio to calculate.
}
\usage{
byGroupHR(survData.list, naive = T, causespec = T, subdist = T, tofile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData.list}{
List of survData arrays
}
  \item{naive}{
Time-independent regression.
}
  \item{causespec}{
Cause-specific Cox regression.
}
  \item{subdist}{
Fine & Grey subdistribution Cox regression.
}
  \item{tofile}{
Should the output be sent to file.
}
}
\details{
The time-independent regression does not account for time of infection (or generally the timing of an intermediate event); the cause-specific Cox regression considers the competing risks to be censored times and the Fine & Grey subdistribution approach retains the competing risk in the risk set.
}
\value{For each element in the list returns

  \item{timeindnaive}{Time-independent model}
  \item{timedeptcausespec}{Time-dependent cause-specific model}
  \item{timedeptsubdistn}{Time-dependent Fine & Grey subdistribution model}
  
Within \code{timeindnaive},

  \item{fit.disch}{Cox regression results with discharge alive event of interest}
  \item{fit.dead}{Cox regression results with death in-hospital event of interest}
  \item{fit.both}{Cox regression results with combined event of interest}
  
Within \code{timedeptcausespec} and \code{timedeptsubdistn}

  \item{fit.alive.timeonly}{Discharge alive as event of interest; only time of infection as covariate}
  \item{fit.alive.full}{Discharge alive as event of interest; time of infection and other risk factors as covariates}
  \item{fit.dead.timeonly}{Death in hospital as event of interest; only time of infection as covariate}
  \item{fit.dead.full}{Death in hospital as event of interest; time of infection and other risk factors as covariates}
  \item{fit.both.timeonly}{Combined event of interest; only time of infection as covariate}
  \item{fit.both.full}{Combined event of interest; time of infection and other risk factors as covariates}
  
Within each of these, \command{coxph} return values.

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
\code{\link{hr.naive}}, \code{\link{hr.tv.causespecific}}, \code{\link{hr.tv.subdistribution}}, \code{\link{coxph}},
\code{\link{byGroupHR.boxplot}}, \code{\link{HRboxplot.batch}}
}
\examples{
require(survival)
data(survDataByGender)
foo <- byGroupHR(survDataByGender)
foo$Male$timedeptcausespec$fit.alive.full


## The function is currently defined as
function (survData.list, naive = T, causespec = T, subdist = T, 
    tofile = FALSE) 
{
    startTime <- Sys.time()
    namesGroup <- names(survData.list)
    results <- list()
    progressbar <- txtProgressBar(min = 0, max = length(namesGroup), 
        style = 3)
    counter = 0
    for (i in namesGroup) {
        counter = counter + 1
        setTxtProgressBar(progressbar, counter)
        timeindnaive <- timedeptcausespec <- timedeptsubdistn <- NA
        if (naive == T) {
            timeindnaive = hr.naive(survData.list[[i]])
        }
        if (causespec == T) {
            timedeptcausespec = hr.tv.causespecific(survData.list[[i]])
        }
        if (subdist == T) {
            timedeptsubdistn = hr.tv.subdistribution(survData.list[[i]])
        }
        results[[i]] <- list(timeindnaive = timeindnaive, timedeptcausespec = timedeptcausespec, 
            timedeptsubdistn = timedeptsubdistn)
    }
    close(progressbar)
    cat(Sys.time() - startTime)
    if (tofile == TRUE) {
        sink(".\\output\\byOrgamismHR.txt")
        results
        sink()
    }
    results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
