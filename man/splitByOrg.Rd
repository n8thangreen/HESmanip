\name{splitByOrg}
\alias{splitByOrg}
\title{Stratify the patient sample by organism type}
\description{
Deprecated by the more general \code{splitByGroup}.
}
\usage{
splitByOrg(total.data, survData.all)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
The complete combined raw data.
}
  \item{survData.all}{
The complete survival data format array with multiple same-patient in-hospital infection records removed.
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
\code{\link{splitByGroup}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (total.data, survData.all) 
{
    survDataByGroup <- list()
    namesGroup <- unique(total.data$lab_group[!is.na(total.data$lab_group)])
    for (name in namesGroup) {
        groupRows <- total.data$lab_group \%in\% c(name, NA)
        uniqueGroupRows <- rmOrganism(total.data[groupRows, ])
        survDataByGroup[[name]] <- survData.all[groupRows, ][uniqueGroupRows, 
            ]
    }
    gram <- c(-1, 1)
    for (i in seq_along(gram)) {
        groupRows <- total.data$Gram \%in\% c(gram[i], NA)
        uniqueGroupRows <- rmOrganism(total.data[groupRows, ])
        survDataByGroup[[as.character(gram[i])]] <- survData.all[groupRows, 
            ][uniqueGroupRows, ]
    }
    survDataByGroup
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
