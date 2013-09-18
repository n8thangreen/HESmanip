\name{summaryTableRF}
\alias{summaryTableRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summaryTableRF}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summaryTableRF(survData)
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
    out <- NULL
    out <- rbind(out, c(NA, NA, NA, NA))
    x <- as.data.frame(with(survData, table(infstatus, cancer)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        cancer)), 1))
    out <- rbind(out, c(x[x$cancer == T, "Freq"], format(y[y$cancer == 
        T, "Freq"], digits = 2)))
    x <- as.data.frame(with(survData, table(infstatus, prem)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        prem)), 1))
    out <- rbind(out, c(x[x$prem == T, "Freq"], format(y[y$prem == 
        T, "Freq"], digits = 2)))
    x <- as.data.frame(with(survData, table(infstatus, cong)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        cong)), 1))
    out <- rbind(out, c(x[x$cong == T, "Freq"], format(y[y$cong == 
        T, "Freq"], digits = 2)))
    x <- as.data.frame(with(survData, table(infstatus, surgical)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        surgical)), 1))
    out <- rbind(out, c(x[x$surgical == T, "Freq"], format(y[y$surgical == 
        T, "Freq"], digits = 2)))
    x <- as.data.frame(with(survData, table(infstatus, cath)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        cath)), 1))
    out <- rbind(out, c(x[x$cath == T, "Freq"], format(y[y$cath == 
        T, "Freq"], digits = 2)))
    x <- as.data.frame(with(survData, table(infstatus, highRisk)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        highRisk)), 1))
    out <- rbind(out, c(x[x$highRisk == T, "Freq"], format(y[y$highRisk == 
        T, "Freq"], digits = 2)))
    out <- rbind(out, c(NA, NA, NA, NA))
    elective <- c("Elective - booked", "Elective - planned", 
        "Elective - from waiting list")
    x <- as.data.frame(with(survData, table(infstatus, hes_admimethdescription \%in\% 
        elective)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        hes_admimethdescription \%in\% elective)), 1))
    out <- rbind(out, c(x[x$Var2 == T, "Freq"], format(y[y$Var2 == 
        T, "Freq"], digits = 2)))
    emergency <- c("Emergency - other means, including patients who arrive via A&E department of another HC provider", 
        "Emergency - via A&E services, including casualty department of provider", 
        "Emergency - via General Practitioner (GP)", "Emergency - via Bed Bureau, including Central Bureau", 
        "Emergency - via consultant out-patient clinic")
    x <- as.data.frame(with(survData, table(infstatus, hes_admimethdescription \%in\% 
        emergency)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        hes_admimethdescription \%in\% emergency)), 1))
    out <- rbind(out, c(x[x$Var2 == T, "Freq"], format(y[y$Var2 == 
        T, "Freq"], digits = 2)))
    out <- rbind(out, aggregate(as.numeric(tail(out, 2)), by = list(c(1, 
        1, 2, 2, 3, 3, 4, 4)), sum)$x)
    neocare <- c("Level 1 intensive care", "Level 2 intensive care")
    x <- as.data.frame(with(survData, table(infstatus, hes_neocaredescription \%in\% 
        neocare)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        hes_neocaredescription \%in\% neocare)), 1))
    out <- rbind(out, c(x[x$Var2 == T, "Freq"], format(y[y$Var2 == 
        T, "Freq"], digits = 2)))
    out <- rbind(out, c(NA, NA, NA, NA))
    x <- as.data.frame(with(survData, table(infstatus, hes_admimethdescription == 
        "Transfer of any admitted patient from another hospital provider")))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        hes_admimethdescription == "Transfer of any admitted patient from another hospital provider")), 
        1))
    out <- rbind(out, c(x[x$Var2 == T, "Freq"], format(y[y$Var2 == 
        T, "Freq"], digits = 2)))
    residence <- c("The usual place of residence, including no fixed abode", 
        "Temporary place of residence when usually resident elsewhere")
    x <- as.data.frame(with(survData, table(infstatus, hes_admisorcdescription \%in\% 
        residence)))
    y <- as.data.frame(prop.table(with(survData, table(infstatus, 
        hes_admisorcdescription \%in\% residence)), 1))
    out <- rbind(out, c(x[x$Var2 == T, "Freq"], format(y[y$Var2 == 
        T, "Freq"], digits = 2)))
    out <- rbind(out, aggregate(as.numeric(tail(out, 2)), by = list(c(1, 
        1, 2, 2, 3, 3, 4, 4)), sum)$x)
    if (ncol(out) == 4) {
        colnames(out) <- c("Non-HA-BSI", "HA-BSI", "\%", "\%")
        out <- out[, c(1, 3, 2, 4)]
    }
    else {
        colnames(out) <- c("Non-HA-BSI", "\%")
    }
    rownames(out) <- c("Risk Factors", "Cancer", "Premature birth", 
        "Congenital disorder", "Surgical", "In-dwelling catheter", 
        "At least one", "Type of Admission", "Elective", "Emergency", 
        "TOTAL", "Intensive neonatal care", "Origin of patient", 
        "Another hospital", "Residence", "TOTAL")
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
