\name{LOSboot.plot}
\alias{LOSboot.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of stay envelope plot
}
\description{
Additionally plots standard error or max min boundaries. Also plots mean or median estimates.
}
\usage{
LOSboot.plot(survData, rlstns = 20, est = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{rlstns}{
Number of bootstrap realisations
}
  \item{est}{Estimator type, either mean or median
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
 \code{\link{LOSboot.plot.wide}}
}
\examples{
data(survDataByGender_sample)
LOSboot.plot(survDataByGender_sample$Male, rlstns=10)


## The function is currently defined as
function (survData, rlstns = 20, est = "mean") 
{
    survData <- survData[!is.na(survData$time) & survData$time > 
        0, ]
    survData <- survData[survData$time >= min(survData$spectime[survData$infstatus == 
        1]), ]
    cases.inf <- survData.list[[group]]$infstatus == 1
    cases.mix <- survData.list[[group]]$infstatus == 0
    ninf <- sum(cases.inf)
    nmix <- sum(cases.mix)
    tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), 
        as.character(0:2)))
    tra[1, 2:3] <- TRUE
    tra[2, 3] <- TRUE
    cLOS <- NULL
    datalist.inf <- datalist.mix <- list()
    for (k in 1:rlstns) {
        sample.mix.new <- survData[cases.mix, ][sample(nmix, 
            replace = T), ]
        sample.inf.new <- survData[cases.inf, ][sample(ninf, 
            replace = T), ]
        survData.new <- rbind(sample.mix.new, sample.inf.new)
        msm.data <- etmFormat(survData.new, type = "")
        if (est == "mean") {
            mvna.data <- mvna::mvna(data = msm.data, state.names = c("0", 
                "1", "2"), tra = tra, cens.name = "cens")
            etm.data <- etm::etm(data = msm.data, state.names = c("0", 
                "1", "2"), tra = tra, cens.name = "cens", s = 0)
            clos.data <- etm::clos(etm.data, aw = TRUE)
            cLOS[k] <- clos.data$e.phi
            datalist.inf[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.case, infstatus = 1)
            datalist.mix[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.control, infstatus = 0)
        }
        if (est == "median") {
            data.wide <- long2wide(msm.data)
            clos.data <- cLOS(my.data = data.wide)
            names(clos.data) <- c("e.phi", "e.phi.med", "times", 
                "phi.case", "phi.control", "phi.case.med", "phi.control.med", 
                "weights", "matrices")
            cLOS[k] <- clos.data$e.phi.med
            datalist.inf[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.case.med, infstatus = 1)
            datalist.mix[[k]] <- cbind(time = clos.data$time, 
                phi = clos.data$phi.control.med, infstatus = 0)
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
    maxtime <- 50
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
        geom_ribbon(alpha = 0.4)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
