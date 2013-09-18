\name{long2wide}
\alias{long2wide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert length of stay data from "new" long format to "old" wide format}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
long2wide(data.long)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.long}{
%%     ~~Describe \code{data.long} here~~
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
#library(changeLOS)
#data(los.data)

#data.long <- prepare.los.data(los.data) ## long format->wide format

data(dat.icu.pneu)
head(dat.icu.pneu)

data.wide <- long2wide(dat.icu.pneu)
head(data.wide)



## The function is currently defined as
function (data.long) 
{
    data.wide <- data.frame(adm.id = sort(unique(data.long$id)))
    data.wide[, c("j.01", "j.02", "j.03", "j.12", "j.13", "cens")] <- Inf
    if ("time" \%in\% names(data.long)) {
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 1], "j.01"] <- data.long$time[data.long$from == 
            0 & data.long$to == 1]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 2], "j.02"] <- data.long$time[data.long$from == 
            0 & data.long$to == 2]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 3], "j.03"] <- data.long$time[data.long$from == 
            0 & data.long$to == 3]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            1 & data.long$to == 2], "j.12"] <- data.long$time[data.long$from == 
            1 & data.long$to == 2]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            1 & data.long$to == 3], "j.13"] <- data.long$time[data.long$from == 
            1 & data.long$to == 3]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$to == 
            "cens"], "cens"] <- data.long$time[data.long$to == 
            "cens"]
    }
    if ("exit" \%in\% names(data.long)) {
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 1], "j.01"] <- data.long$exit[data.long$from == 
            0 & data.long$to == 1]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 2], "j.02"] <- data.long$exit[data.long$from == 
            0 & data.long$to == 2]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            0 & data.long$to == 3], "j.03"] <- data.long$exit[data.long$from == 
            0 & data.long$to == 3]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            1 & data.long$to == 2], "j.12"] <- data.long$exit[data.long$from == 
            1 & data.long$to == 2]
        data.wide[data.wide$adm.id \%in\% data.long$id[data.long$from == 
            1 & data.long$to == 3], "j.13"] <- data.long$exit[data.long$from == 
            1 & data.long$to == 3]
        for (i in data.long$id[data.long$to == "cens"]) {
            if (sum(data.long$id == i) > 1) {
                data.wide[data.wide$adm.id == i, "cens"] <- data.long$exit[data.long$id == 
                  i][2]
            }
            else {
                data.wide[data.wide$adm.id == i, "cens"] <- data.long$exit[data.long$id == 
                  i][1]
            }
        }
    }
    data.wide[data.wide$j.12 == Inf & data.wide$j.02 == Inf & 
        data.wide$j.03 == Inf & data.wide$cens == Inf & data.wide$j.13 == 
        Inf & data.wide$j.01 != Inf, "cens"] <- data.wide[data.wide$j.12 == 
        Inf & data.wide$j.02 == Inf & data.wide$j.03 == Inf & 
        data.wide$cens == Inf & data.wide$j.13 == Inf & data.wide$j.01 != 
        Inf, "j.01"]
    data.wide
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
