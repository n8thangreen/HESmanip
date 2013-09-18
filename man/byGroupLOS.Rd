\name{byGroupLOS}
\alias{byGroupLOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of stay by groups}
\description{
Batch calculate all the cumulative hazards, transition probabilites and length of stay for each of the organism groups in turn
}
\usage{
byGroupLOS(survData.list, type, standerr = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData.list}{
%%     ~~Describe \code{survData.list} here~~
}
  \item{type}{
Model type: cause-specific "", subdistribution with discharge alive as event of interest "alive" or subdistribution with in-hospital death as the event of interest "dead".
}
  \item{standerr}{
Should the standard error be calculated. This may take some minutes.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns list of lists. Below initial sample stratification:

  \item{mvna.data}{Output from \command{mvna} contained in \pkg{mvna}}
  \item{clos.data}{Excess length of stay output}
  \item{se}{Bootstrap standard error estimate}
  \item{etm.data}{Output from \command{etm} contained in \pkg{etm}}
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
\code{\link{byGroupLOS.boxplot}}
}
\examples{
require(mvna)
require(etm)
data(survDataByGender)

foo <- byGroupLOS(survDataByGender,"")
foo$Male$clos.data



## The function is currently defined as
function (survData.list, type, standerr = FALSE) 
{
    startTime <- Sys.time()
    namesGroup <- names(survData.list)
    progressbar <- txtProgressBar(min = 0, max = length(namesGroup), 
        style = 3)
    counter = 0
    tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), 
        as.character(0:2)))
    tra[1, 2:3] <- TRUE
    tra[2, 3] <- TRUE
    results <- list()
    for (i in namesGroup) {
        counter = counter + 1
        setTxtProgressBar(progressbar, counter)
        if (i == "") {
            stop("no group name")
        }
        if (nrow(survData.list[[i]]) == 0) {
            results[[i]] <- list(NA, NA, NA)
        }
        else {
            results[[i]] <- listLOS(survData.list[[i]], tra, 
                type, standerr)
        }
    }
    close(progressbar)
    cat(Sys.time() - startTime)
    results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
