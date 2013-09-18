\name{artifCens}
\alias{artifCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Artificial censoring}
\description{
Artificially impose censoring at either earliest censoring time or a user defined time.
}
\usage{
artifCens(total.data, censtime = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
%%     ~~Describe \code{total.data} here~~
}
  \item{censtime}{
Fraction of time between smallest and largest administrative censoring times [0,1]
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
%%  ~~who you are~~
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
function (total.data, censtime = NA) 
{
    mix.disdteMax <- max(total.data$hes_disdte[total.data$infstatus == 
        0], na.rm = TRUE)
    inf.disdteMax <- max(total.data$hes_disdte[total.data$infstatus == 
        1], na.rm = TRUE)
    artifCens.time <- min(mix.disdteMax, inf.disdteMax)
    if (is.numeric(censtime) & !is.na(censtime)) {
        artifCens.time <- artifCens.time + (max(mix.disdteMax, 
            inf.disdteMax) - artifCens.time) * censtime
    }
    artifCens.dis <- total.data$hes_disdte > artifCens.time
    total.data$hes_dismethdescription[artifCens.dis] <- "Not applicable patient still in hospital"
    total.data$hes_dischargeddeadoralive[artifCens.dis] <- ""
    total.data$hes_disdte[artifCens.dis] <- artifCens.time
    total.data$missingType[artifCens.dis] <- 0
    artifCens.spec <- total.data$lab_Specimendate > artifCens.time
    total.data$infstatus[artifCens.spec] <- 0
    total.data$adm_diff_specdate[artifCens.spec] <- total.data$dis_diff_specdate[artifCens.spec] <- total.data$lab_Specimendate[artifCens.spec] <- NA
    total.data$adm_diff_dis <- with(total.data, hes_disdte - 
        hes_admdte)
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
