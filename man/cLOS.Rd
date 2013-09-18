\name{cLOS}
\alias{cLOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extension of clos() with medians}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cLOS(x = NA, my.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{my.data}{Wide format j.01, j.02, j.03, j.12, j.13, cens}
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
\code{\link{clos}}
}
\examples{
library(etm)
library(bootstrap)
library(changeLOS)
data(los.data)

data.long <- prepare.los.data(los.data) ## long format->wide format

tra <- matrix(TRUE, 4, 4, dimnames = list(as.character(0:3), as.character(0:3)))
diag(tra)<-FALSE
tra[2:4,1] <- tra[3:4,2] <- tra[4,3] <- tra[3,4] <- FALSE

etm.data <- etm(data.long, c("0","1","2","3"), tra, "cens", s=0)
summary(etm.data)
xyplot(etm.data)


clos.long <- clos(etm.data, aw=TRUE)       # mean length of stay

clos.wide <- cLOS(my.data=los.data)       # median (and mean) length of stay






## The function is currently defined as
function (x = NA, my.data) 
{
    require(survival)
    library(bootstrap)
    if (is.na(x)) {
        x <- 1:length(my.data[, 1])
    }
    my.data <- my.data[x, ]
    my.data$cens.0 <- my.data$cens
    my.data$cens.0[is.finite(my.data$j.01)] <- Inf
    my.data$cens.1 <- my.data$cens
    my.data$cens.1[is.infinite(my.data$j.01)] <- Inf
    jump.times <- sort(unique(c(my.data$j.01, my.data$j.02, my.data$j.03, 
        my.data$j.12, my.data$j.13)))
    jump.times <- jump.times[is.finite(jump.times)]
    jump.matrices <- array(0, c(4, 4, length(jump.times)))
    for (i in 1:length(jump.times)) {
        jump.matrices[1, 2, i] <- length(my.data$j.01[my.data$j.01 == 
            jump.times[i]])
        jump.matrices[1, 3, i] <- length(my.data$j.02[my.data$j.02 == 
            jump.times[i]])
        jump.matrices[1, 4, i] <- length(my.data$j.03[my.data$j.03 == 
            jump.times[i]])
        jump.matrices[2, 3, i] <- length(my.data$j.12[my.data$j.12 == 
            jump.times[i]])
        jump.matrices[2, 4, i] <- length(my.data$j.13[my.data$j.13 == 
            jump.times[i]])
        risk.0 <- length(my.data[, 1]) - length(c(my.data$j.01, 
            my.data$j.02, my.data$j.03, my.data$cens.0)[c(my.data$j.01, 
            my.data$j.02, my.data$j.03, my.data$cens.0) < jump.times[i]])
        risk.1 <- length(my.data$j.01[my.data$j.01 < jump.times[i]]) - 
            length(c(my.data$j.12, my.data$j.13, my.data$cens.1)[c(my.data$j.12, 
                my.data$j.13, my.data$cens.1) < jump.times[i]])
        if (risk.0 > 0) 
            jump.matrices[1, , i] <- jump.matrices[1, , i]/risk.0
        if (risk.1 > 0) 
            jump.matrices[2, , i] <- jump.matrices[2, , i]/risk.1
        jump.matrices[, , i] <- jump.matrices[, , i] + diag(1, 
            4, 4) - diag(apply(jump.matrices[, , i], 1, sum))
    }
    my.times <- sort(unique(c(jump.times, max(my.data$cens[is.finite(my.data$cens)], 
        jump.times))))
    los <- matrix(data = rep(my.times, 3), ncol = 3, byrow = FALSE, 
        dimnames = list(NULL, c("Time", "Given in state 1", "Given in state 0")))
    los[length(my.times) - 1, 2:3] <- rep(max(my.times), 2)
    med.los <- matrix(data = rep(my.times, 3), ncol = 3, byrow = FALSE, 
        dimnames = list(NULL, c("Time", "Given in state 1", "Given in state 0")))
    aj <<- array(NA, c(4, 4, 1))
    aj[, , 1] <- diag(1, 4, 4)
    my.function <- function(x, y) {
        x \%*\% aj[, , y]
    }
    for (i in (length(my.times) - 2):1) {
        diffs <- diff(my.times[(i + 1):length(my.times)])
        my.matrix <- jump.matrices[, , length(jump.times[jump.times <= 
            my.times[i + 1]])]
        aj <- array(apply(X = diag(1:dim(aj)[3]), 1, my.function, 
            x = my.matrix), c(4, 4, dim(aj)[3]))
        los[, 2][i] <- my.times[i + 1] + matrix(diffs, nrow = 1) \%*\% 
            matrix(aj[2, 2, ], ncol = 1)
        los[, 3][i] <- my.times[i + 1] + matrix(diffs, nrow = 1) \%*\% 
            matrix((aj[1, 1, ] + aj[1, 2, ]), ncol = 1)
        ecdf.inf <- aj[2, 3, ] + aj[2, 4, ]
        med.min <- min(which(ecdf.inf >= 0.5))
        med.los[, 2][i] <- my.times[i + med.min]
        ecdf.adm <- aj[1, 3, ] + aj[1, 4, ]
        med.min <- min(which(ecdf.adm >= 0.5))
        med.los[, 3][i] <- my.times[i + med.min]
        aj <- array(c(diag(1, 4, 4), aj), c(4, 4, (dim(aj)[3] + 
            1)))
    }
    T0 <- Surv(apply(as.matrix(my.data[, c("j.01", "j.02", "j.03", 
        "cens")]), 1, min), 1 - is.finite(my.data$cens.0))
    T0.fit <- survfit(T0 ~ 1)
    T0.fit$time <- T0.fit$time[T0.fit$n.event != 0]
    T0.fit$surv <- T0.fit$surv[T0.fit$n.event != 0]
    my.weights <- diff(c(0, 1 - T0.fit$surv))
    estimate <- matrix((los[, 2] - los[, 3])[is.element(los[, 
        1], T0.fit$time)], nrow = 1) \%*\% matrix(my.weights, ncol = 1)
    diff.e <- matrix((med.los[, 2] - med.los[, 3])[is.element(med.los[, 
        1], T0.fit$time)], nrow = 1)
    nNA <- !is.na(diff.e)
    diff.e <- diff.e[nNA]
    estimate.med <- diff.e \%*\% matrix(my.weights[nNA], ncol = 1)
    my.return <- list(cLOS = estimate, cLOS.med = estimate.med, 
        times = jump.times, e.given.1 = c(los[, 2]), e.given.0 = c(los[, 
            3]), e.given.1.med = c(med.los[, 2]), e.given.0.med = c(med.los[, 
            3]), weights = my.weights, matrices = jump.matrices)
    return(my.return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
