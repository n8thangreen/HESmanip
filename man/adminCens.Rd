\name{adminCens}
\alias{adminCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Administrative censoring}
\description{
If administratively censored fill-in censoring time and add flag. 
}
\usage{
adminCens(total.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.data}{
%%     ~~Describe \code{total.data} here~~
}
}
\details{
Infected and non-infected datasets handled separately. Maximum discharge/death date for each sample separately i.e. administrative censoring times (not necessarily equal).
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
    mix.disdteMax <- max(total.data$hes_disdte[total.data$infstatus == 
        0], na.rm = TRUE)
    inf.disdteMax <- max(total.data$hes_disdte[total.data$infstatus == 
        1], na.rm = TRUE)
    inhosp <- total.data$missingTime == 0 & total.data$hes_dismethdescription \%in\% 
        c("Not applicable patient still in hospital", "Not known")
    total.data$hes_disdte[total.data$infstatus == 0 & inhosp] <- mix.disdteMax
    total.data$hes_disdte[total.data$infstatus == 1 & inhosp] <- inf.disdteMax
    total.data$missingTime <- total.data$missingTime | inhosp
    total.data$missingType[inhosp] <- 0
    total.data$adm_diff_dis <- ifelse(total.data$hes_disdte - 
        total.data$hes_admdte == 0, 0.5, total.data$hes_disdte - 
        total.data$hes_admdte)
    total.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
