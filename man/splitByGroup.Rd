\name{splitByGroup}
\alias{splitByGroup}
\title{Split patient record data in subsets which all contain all non-infected patient records}
\description{
This is a generalisation of the \code{splitByOrg} function and has in effect replaced it.
}
\usage{
splitByGroup(total.data, survData.total, groupCol = "lab_group")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
The complete combined raw data.
}
  \item{survData.total}{
The complete survival data format array with multiple same-patient in-hospital infection records retained.
}
  \item{groupCol}{The field name to split into lists by. Currently use \code{lab_group} for organism stratification, \code{Gram}, \code{CombinedMET} or \code{VAN}.  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of stratified survival data patient records arrays.}
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
\code{\link{splitByOrg}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (total.data, survData.total, groupCol = "lab_group") 
{
    survDataByGroup <- list()
    namesGroup <- as.character(unique(total.data[!is.na(total.data[, 
        groupCol]), groupCol]))
    for (name in namesGroup) {
        groupRows <- (total.data[, groupCol] == name & total.data$infstatus == 
            1) | total.data$infstatus == 0
        groupRows[is.na(groupRows)] <- FALSE
        uniqueGroupRows <- rmOrganism(total.data[groupRows, ])
        survDataByGroup[[name]] <- survData.total[groupRows, 
            ][uniqueGroupRows, ]
    }
    survDataByGroup
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
