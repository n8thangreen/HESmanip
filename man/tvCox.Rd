\name{tvCox}
\alias{tvCox}
\title{
Time-varying regression data representation
}
\description{
Time-dependent hazard long array representation of time-to-event data.
}
\usage{
tvCox(survData, covariates, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
%%     ~~Describe \code{survData} here~~
}
  \item{covariates}{
vector of covariate names (column headings) in \code{survData}.
}
  \item{type}{
which model to use, cause-specific or subdistribution (Fine & Gray)
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
function (survData, covariates, type) 
{
    if (type == "alive") {
        time <- survData$time2disch
    }
    else if (type == "death") {
        time <- survData$time2death
    }
    else if (type == "") {
        time <- survData$time
    }
    rows.inf <- survData$infstatus == 1
    survData.inf <- survData[rows.inf, ]
    survData.mix <- survData[!rows.inf, ]
    ninf <- rep(1:nrow(survData.inf), each = 2)
    epsilon <- 0.5
    time.inf <- time[rows.inf]
    time.mix <- time[!rows.inf]
    originalid <- 0
    tvdata.inf <- data.frame(id = ninf, tstart = NA, tstop = NA, 
        inf = NA, disch = NA, death = NA, survData.inf[ninf, 
            c("infstatus", covariates)])
    numRows.new <- nrow(tvdata.inf)
    if (numRows.new > 0) {
        for (i in seq(1, numRows.new, by = 2)) {
            originalid <- originalid + 1
            spectime <- survData.inf$spectime[originalid]
            time.inf.id <- time.inf[originalid]
            tvdata.inf[i, "tstart"] <- 0
            tvdata.inf[i, "tstop"] <- ifelse(spectime == 0, epsilon, 
                spectime)
            tvdata.inf[i, "disch"] <- 0
            tvdata.inf[i, "inf"] <- 0
            tvdata.inf[i, "death"] <- 0
            tvdata.inf[i + 1, "tstart"] <- tvdata.inf[i, "tstop"]
            tvdata.inf[i + 1, "tstop"] <- ifelse(time.inf.id == 
                spectime, time.inf.id + epsilon, time.inf.id)
            tvdata.inf[i + 1, "disch"] <- (1 - as.numeric(survData.inf$event[originalid])) & 
                survData.inf$missingType[originalid]
            tvdata.inf[i + 1, "inf"] <- 1
            tvdata.inf[i + 1, "death"] <- unlist(survData.inf$event[originalid]) & 
                survData.inf$missingType[originalid]
        }
    }
    tvdata.mix <- data.frame(id = originalid + 10 + seq_along(survData.mix$time), 
        tstart = 0, tstop = ifelse(time.mix > 0, time.mix, epsilon), 
        inf = 0, disch = (1 - as.numeric(survData.mix[, "event"])) & 
            survData.mix[, "missingType"], death = unlist(survData.mix[, 
            "event"]) & survData.mix[, "missingType"], survData.mix[, 
            c("infstatus", covariates)])
    tvdata <- rbind(tvdata.inf, tvdata.mix)
    tvdata <- cbind(tvdata, status = tvdata$disch + (2 * tvdata$death))
    tvdata <- tvdata[!is.na(tvdata[, "tstop"]), ]
    names(tvdata)[which(names(tvdata) == "gender")] <- "sex"
    names(tvdata)[which(names(tvdata) == "infstatus")] <- "originf"
    tvdata
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
