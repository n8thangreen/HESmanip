\name{getSurvData}
\alias{getSurvData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform to survival data format}
\description{
Create a simpler survival analysis format array; basically a subset of the total.data columns with the dates removed so all time origins are equal.
}
\usage{
getSurvData(total.data)
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
function (total.data) 
{
    survData <- within(total.data, {
        time <- as.numeric(adm_diff_dis)
        event <- hes_dischargeddeadoralive == "Death"
        age <- hes_ageatstartofepisode
        gender <- droplevels(hes_sexdescription)
        spectime <- as.numeric(adm_diff_specdate)
        agegr <- cut(age, breaks = c(0, 2, 5, 10, 18, 40, 60, 
            100, 10000))
        rm(lab_opieid, hes_ID, lab_Specimendate, hes_admdte, 
            hes_disdte, adm_diff_specdate, dis_diff_specdate, 
            hes_dismethdescription, adm_diff_dis, hes_sexdescription, 
            hes_ageatstartofepisode, hes_dischargeddeadoralive, 
            lab_organismname, lab_group, Gram)
    })
    survData$age[survData$age > 1000] <- 0
    survData$agegr <- factor(survData$agegr, c("(100,1e+04]", 
        "(0,2]", "(2,5]", "(5,10]", "(10,18]", "(18,40]", "(40,60]", 
        "(60,100]"))
    survData <- subdistn(survData)
    survData
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
