\name{LOStable}
\alias{LOStable}
\title{Hospital length of stay table}
\description{
Stratified data estimates of length of stay with standard errors.
}
\usage{
LOStable(output.LOS, se = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.LOS}{
%%     ~~Describe \code{output.LOS} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
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
\code{\link{listLOS}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (output.LOS, se = TRUE) 
{
    LOS <- NULL
    dp <- 2
    for (i in names(output.LOS)) {
        if (is.na(output.LOS[[i]][[1]][1])) {
            LOS <- c(LOS, NA)
        }
        else {
            if (se) {
                LOS <- c(LOS, paste(round(output.LOS[[i]]$clos.data$e.phi, 
                  dp), " (", round(output.LOS[[i]]$se, dp), ")", 
                  sep = ""))
            }
            else {
                LOS <- c(LOS, paste(round(output.LOS[[i]]$clos.data$e.phi, 
                  dp), " [", round(output.LOS[[i]]$clos.data$e.phi - 
                  1.96 * output.LOS[[i]]$se, dp), ", ", round(output.LOS[[i]]$clos.data$e.phi + 
                  1.96 * output.LOS[[i]]$se, dp), "]", sep = ""))
            }
        }
    }
    res <- data.frame(Group = names(output.LOS), `Excess LOS` = LOS)
    rownames(res) <- names(output.LOS)
    res <- res[, -1, drop = FALSE]
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
