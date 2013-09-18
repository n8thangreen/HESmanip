\name{findDuplicates}
\alias{findDuplicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find patient record duplicates}
\description{
Find (and remove) duplicate same-patient entries on same-day of admission with missing event-type and observed (dead, alive or administrative censoring) because these have the same admission and discharge data as one another; making sure to retain all relevant information.
}
\usage{
findDuplicates(total.data, nCodesInd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{Full combined infected and non-infected patients data set}
  \item{nCodesInd}{The names of the risk factor or comorbidity column names in total.dat}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (total.data, nCodesInd) 
{
    uniqueid <- unique(total.data$hes_ID)
    progressbar <- txtProgressBar(min = 0, max = length(uniqueid), 
        style = 3)
    counter = 0
    for (i in uniqueid) {
        counter = counter + 1
        setTxtProgressBar(progressbar, counter)
        ind <- which(total.data$hes_ID == i)
        uniquetime <- unique(total.data$hes_admdte[ind])
        for (j in uniquetime) {
            subrows <- ind[total.data$hes_admdte[ind] == j]
            if (length(subrows) > 1) {
                total.data[subrows, nCodesInd] <- matrix(apply(total.data[subrows, 
                  nCodesInd], MARGIN = 2, FUN = any), nrow = length(subrows), 
                  ncol = length(nCodesInd), byrow = TRUE)
                if (any(!is.na(total.data$lab_Specimendate[subrows]))) {
                  total.data[subrows, ] <- total.data[subrows, 
                    ][order(total.data$lab_Specimendate[subrows]), 
                    ]
                }
                if (any(is.na(total.data[subrows, "lab_opieid"])) | 
                  anyDuplicated(total.data[subrows, "lab_opieid"])) {
                  total.data <- rmDuplicates(total.data, subrows)
                }
            }
        }
    }
    close(progressbar)
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
