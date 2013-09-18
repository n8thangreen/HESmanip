\name{dat.icu.pneu}
\alias{dat.icu.pneu}
\docType{data}
\title{Sample ICU pneumonia time to event data from kmi package in "long" format}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(dat.icu.pneu)}
\format{
  A data frame with 1421 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{entry}}{a numeric vector}
    \item{\code{exit}}{a numeric vector}
    \item{\code{from}}{a factor with levels \code{0} \code{1}}
    \item{\code{to}}{a factor with levels \code{1} \code{2} \code{cens}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(dat.icu.pneu)
## maybe str(dat.icu.pneu) ; plot(dat.icu.pneu) ...
}
\keyword{datasets}
