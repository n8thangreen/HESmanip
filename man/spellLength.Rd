\name{spellLength}
\alias{spellLength}
\title{Find actual start and end of a patients in-hospital spell}
\description{

}
\usage{
spellLength(adm, disch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adm}{
Admission dates for (subset of) patient records.
}
  \item{disch}{
Discharge dates for (subset of) patient record.
}
}
\details{
Chains of linked discharge and admission dates are made to then determine what the overall first admission and last disscharge times are and then these are substituted in to all of the relevant records.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (adm, disch) 
{
    for (i in 1:length(adm)) {
        dupl <- adm[i] == disch
        disch[dupl] <- disch[i]
        if (any(dupl, na.rm = TRUE)) {
            adm[i] <- adm[dupl] <- min(adm[dupl], na.rm = TRUE)
        }
    }
    list(adm = adm, disch = disch)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
