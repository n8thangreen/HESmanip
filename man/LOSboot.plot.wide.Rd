\name{LOSboot.plot.wide}
\alias{LOSboot.plot.wide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of stay bootstrap plot with envelope bounds for "old" wide event data format}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LOSboot.plot.wide(data.wide, rlstns = 20, est = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.wide}{
wide format of time to event data
}
  \item{rlstns}{
number of bootstrap realisations
}
  \item{est}{
estimator type, either mean or median
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
\code{\link{LOSboot.plot}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data.wide, rlstns = 20, est = "mean") 
{
    cases.inf <- data.wide$j.01 != Inf
    cases.mix <- data.wide$j.01 == Inf
    ninf <- sum(cases.inf)
    nmix <- sum(cases.mix)
    tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), 
        as.character(0:2)))
    tra[1, 2:3] <- TRUE
    tra[2, 3] <- TRUE
    cLOS <- NULL
    datalist.inf <- datalist.mix <- list()
    clos.data.orig <- cLOS(my.data = data.wide)
    for (k in 1:rlstns) {
        sample.mix.new <- data.wide[cases.mix, ][sample(nmix, 
            replace = T), ]
        sample.inf.new <- data.wide[cases.inf, ][sample(ninf, 
            replace = T), ]
        data.wide.new <- rbind(sample.mix.new, sample.inf.new)
        clos.data <- cLOS(my.data = data.wide.new)
        names(clos.data) <- c("e.phi", "e.phi.med", "times", 
            "phi.case", "phi.control", "phi.case.med", "phi.control.med", 
            "weights", "matrices")
        if (est == "median") {
            cLOS[k] <- clos.data$e.phi.med
            datalist.inf[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.case.med, infstatus = 1)
            datalist.mix[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.control.med, infstatus = 0)
        }
        if (est == "mean") {
            cLOS[k] <- clos.data$e.phi
            datalist.inf[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.case, infstatus = 1)
            datalist.mix[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.control, infstatus = 0)
        }
    }
    data.inf <- Reduce(function(x, y) merge(x, y, by = c("time", 
        "infstatus")), datalist.inf)
    data.mix <- Reduce(function(x, y) merge(x, y, by = c("time", 
        "infstatus")), datalist.mix)
    data.inf <- data.inf[order(data.inf$time), ]
    data.mix <- data.mix[order(data.mix$time), ]
    data <- rbind(data.inf, data.mix)
    cols <- 2 + (1:rlstns)
    maxtime <- 500
    data.all <- cbind(data[data$time < maxtime, c(1, 2)], LOS = apply(data[data$time < 
        maxtime, cols], 1, mean), LOS.lo = apply(data[data$time < 
        maxtime, cols], 1, mean) - 1.96 * apply(data[data$time < 
        maxtime, cols], 1, sd), LOS.hi = apply(data[data$time < 
        maxtime, cols], 1, mean) + 1.96 * apply(data[data$time < 
        maxtime, cols], 1, sd))
    data.all <- cbind(data.all, Status = as.factor(data.all$infstatus))
    x11()
    ggplot(data = data.all, aes(x = time, y = LOS, ymin = LOS.lo, 
        ymax = LOS.hi, fill = Status, linetype = Status)) + geom_line() + 
        geom_ribbon(alpha = 0.4) + ylim(0, 75)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
