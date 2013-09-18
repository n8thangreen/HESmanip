\name{summaryTableGroup}
\alias{summaryTableGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summaryTableGroup}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summaryTableGroup(survDataByGroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survDataByGroup}{
%%     ~~Describe \code{survDataByGroup} here~~
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
function (survDataByGroup) 
{
    out <- NA
    for (group in names(survDataByGroup)) {
        out <- rbind(out, c(group, paste(round(mean(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T)), "/", median(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T), sep = ""), round(table(survDataByGroup[[group]]$gender[survDataByGroup[[group]]$infstatus == 
            1])[1]), paste(round(mean(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T)), "/", median(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T), " (", round(sd(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T)), ")", sep = ""), paste(round(mean(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T)), "/", median(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T), " (", round(sd(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus == 
            1], na.rm = T)), ")", sep = ""), round(table(survDataByGroup[[group]]$event[survDataByGroup[[group]]$infstatus == 
            1])[2])))
    }
    colnames(out) <- c("Group", "Age (years)", "Sex", "LOS (days)", 
        "Inf time (days)", "Deaths")
    rownames(out) <- out[, "Group"]
    out <- out[!is.na(out[, 1]), -1]
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
