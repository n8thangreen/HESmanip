\name{summaryTableAll}
\alias{summaryTableAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summaryTableAll}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summaryTableAll(survData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
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
function (survData) 
{
    survData.mix <- survData[survData$infstatus == 0, ]
    survData.inf <- survData[survData$infstatus == 1, ]
    out <- rbind(c(format(nrow(survData.inf), nsmall = 1), round(nrow(survData.inf)/nrow(survData), 
        3), round(nrow(survData.mix)), round(nrow(survData.mix)/nrow(survData), 
        3), round(nrow(survData)), 1), c(paste(round(mean(survData.inf$age)), 
        "/", median(survData.inf$age)), round(sd(survData.inf$age)), 
        paste(round(mean(survData.mix$age)), "/", median(survData.mix$age)), 
        round(sd(survData.mix$age)), paste(round(mean(survData$age)), 
            "/", median(survData$age)), round(sd(survData$age))), 
        c(paste(round(mean(survData.inf$time)), "/", median(survData.inf$time)), 
            round(sd(survData.inf$time)), paste(round(mean(survData.mix$time, 
                na.rm = T)), "/", median(survData.mix$time, na.rm = T)), 
            round(sd(survData.mix$time, na.rm = T)), paste(round(mean(survData$time, 
                na.rm = T)), "/", median(survData$time, na.rm = T)), 
            round(sd(survData$time, na.rm = T))), c(paste(round(mean(survData.inf$spectime)), 
            "/", median(survData.inf$spectime)), round(sd(survData.inf$spectime)), 
            paste(round(mean(survData.mix$spectime)), "/", median(survData.mix$spectime)), 
            round(sd(survData.mix$spectime)), paste(round(mean(survData$spectime)), 
                "/", median(survData$spectime)), round(sd(survData$spectime))), 
        c(round(table(survData.inf$event)[2]), round(table(survData.inf$event)[2]/(table(survData.inf$event)[2] + 
            table(survData.inf$event)[1]), 2), round(table(survData.mix$event)[2]), 
            round(table(survData.mix$event)[2]/(table(survData.mix$event)[2] + 
                table(survData.mix$event)[1]), 2), round(table(survData$event)[2]), 
            round(table(survData$event)[2]/(table(survData$event)[2] + 
                table(survData$event)[1]), 3)), c(round(table(survData.inf$gender)[1]), 
            round(table(survData.inf$gender)[1]/(table(survData.inf$gender)[1] + 
                table(survData.inf$gender)[2]), 2), round(table(survData.mix$gender)[1]), 
            round(table(survData.mix$gender)[1]/(table(survData.mix$gender)[1] + 
                table(survData.mix$gender)[2]), 2), round(table(survData$gender)[1]), 
            round(table(survData$gender)[1]/(table(survData$gender)[1] + 
                table(survData$gender)[2]), 2)))
    rownames(out) <- c("Sample size", "Ages (years)", "LOS (days)", 
        "Inf times (days)", "Deaths", "Sex (F)")
    colnames(out) <- c("HA-BSI", "\%", "Non-HA-BSI", "\%", "All", 
        "\%")
    return(pandoc.table(out, caption = "caption:...", style = "grid", 
        split.tables = Inf, justify = "left"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
