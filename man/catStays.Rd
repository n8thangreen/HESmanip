\name{catStays}
\alias{catStays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Concatenate hospital stays}
\description{
Concatenate same spell records when the patient is moved within hospital i.e. discharge=arrival time
}
\usage{
catStays(total.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
%%     ~~Describe \code{total.data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
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
function (total.data) 
{
    uniqueid <- unique(total.data$hes_ID)
    progressbar <- txtProgressBar(min = 0, max = length(uniqueid), 
        style = 3)
    counter = 0
    for (i in uniqueid) {
        counter = counter + 1
        setTxtProgressBar(progressbar, counter)
        ind <- which(total.data$hes_ID == i)
        newTimes <- spellLength2(total.data$hes_admdte[ind], 
            total.data$hes_disdte[ind])
        total.data$hes_admdte[ind] <- newTimes$adm
        total.data$hes_disdte[ind] <- newTimes$disch
        total.data$missingTime[ind][is.na(total.data$hes_disdte[ind])] <- 0
    }
    close(progressbar)
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
