\name{rmOrganism}
\alias{rmOrganism}
\title{Remove duplicate entries which only differ by HCAI organism type}
\description{
Within an organism group take the first incidence of infection i.e earliest specimen time and remove the rest.
}
\usage{
rmOrganism(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Patient record data with approriate fields.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data) 
{
    data.drop <- subset(data, select = -c(adm_diff_specdate, 
        dis_diff_specdate, lab_Specimendate, lab_opieid, lab_organismname, 
        lab_group, Gram))
    !duplicated(data.drop)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
