\name{antibioticLookup.dat}
\alias{antibioticLookup.dat}
\docType{data}
\title{Antibiotic look-up table}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(antibioticLookup.dat)}
\format{
  A data frame with 59 observations on the following 5 variables.
  \describe{
    \item{\code{source}}{a factor with levels \code{} \code{BNF-C} \code{HPA HCAI} \code{Lambert2011}}
    \item{\code{suspected_organism}}{a character vector}
    \item{\code{class_antibiotic}}{a character vector}
    \item{\code{antibiotic}}{a character vector}
    \item{\code{drug_code}}{a factor with levels \code{} \code{AMI} \code{AMOX} \code{AMP} \code{AMPAMOX} \code{BPEN} \code{CFOX} \code{CHLOR} \code{CLF} \code{CTX} \code{CTZ} \code{FLU} \code{FOR} \code{GENT} \code{IMI} \code{LAT} \code{MER} \code{MERO} \code{MOX} \code{MTR} \code{NET} \code{OXA} \code{PEN} \code{PIPTAZ} \code{STR} \code{TICCLAV} \code{TOB} \code{VAN} \code{van2}}
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
data(antibioticLookup.dat)
## maybe str(antibioticLookup.dat) ; plot(antibioticLookup.dat) ...
}
\keyword{datasets}
