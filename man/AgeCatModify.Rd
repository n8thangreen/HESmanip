\name{AgeCatModify}
\alias{AgeCatModify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Age category modification}
\description{
Find which corresponding age category in the guidelines are the patients in the LabBase2 database.
}
\usage{
AgeCatModify(inf.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf.data}{Infected patient data.}
}
\details{
Find which corresponding age category in the guidelines are the patients in the LabBase2 database.
The age categories different in guidelines and databases so need to regroup and relabel.
  
   7001: 48 hours old;                               
   7002-7003: 48 hours - 1 month old;              
   7004-7007 and 1-18: 1 month - 18 years old;   
   19-200: Over 18 years old                          
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
function (inf.data) 
{
    patient.age_cat <<- cut(inf.data[, "hes_ageatstartofepisode"], 
        breaks = c(0, 18, 200, 7001.5, 7003.5, 7008))
    patient.age_cat <- factor(patient.age_cat, c("(200,7002]", 
        "(7002,7004]", "(7004,7008]", "(0,18]", "(18,200]"))
    patient.age_cat[patient.age_cat == "(0,18]"] <- "(7004,7008]"
    patient.age_cat <- droplevels(patient.age_cat)
    levels(patient.age_cat) <- 1:4
    patient.age_cat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
