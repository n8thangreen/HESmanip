\name{whichPatientDrug.joins}
\alias{whichPatientDrug.joins}
\title{
Identify which patients have anti-microbial resistance and which are sensitive.
}
\description{
Patients are grouped into either tested to be resistant, susceptible or not tested against a specific antibiotic.
}
\usage{
whichPatientDrug.joins(inf.drug.data, drugGuideline.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf.drug.data}{
data for each patient
}
  \item{drugGuideline.data}{
look-up table of guidline drugs for each organism
}
}
\details{
This is an alternative to the original \code{whichPatientDrug} function which uses loops rather than the joining operations.
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

\seealso{
\code{\link{whichPatientDrug}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (inf.drug.data, drugGuideline.data) 
{
    lookup <- drugGuideline.data[, c("suspected_organism", "Generic.Antimicrobial.Code")]
    colnames(lookup) <- c("organism", "drug")
    lookup <- lookup[!duplicated(lookup), ]
    colnames(inf.drug.data) <- drug.names
    inf.drug.data$patient <- rownames(inf.drug.data)
    inf.drug.data <- cbind(inf.drug.data, organism = inf.data$lab_group)
    moltendata <- melt(inf.drug.data, id = c("patient", "organism"))
    colnames(moltendata) <- c("patient", "organism", "drug", 
        "value")
    moltendata <- moltendata[order(moltendata$patient), ]
    lookup <- apply(lookup, 2, toupper)
    moltendata <- apply(moltendata, 2, toupper)
    GLoutput <- merge(moltendata, lookup, by = c("drug", "organism"))
    GLoutput <- subset(GLoutput, select = c("drug", "patient", 
        "value"))
    GLoutput.generic <- merge(moltendata, lookup[lookup[, "organism"] == 
        "", ], by = "drug")
    GLoutput.generic <- subset(GLoutput.generic, select = c("drug", 
        "patient", "value"))
    GLoutput <- rbind(GLoutput, GLoutput.generic)
    GLoutput$patient <- as.numeric(as.character(GLoutput$patient))
    GLoutput <- GLoutput[order(GLoutput$patient), ]
    GLoutput <- GLoutput[!duplicated(GLoutput), ]
    GLstrata <- list()
    GLstrata$res <- GLoutput[GLoutput$value == 1 & !is.na(GLoutput$value), 
        ]
    GLstrata$sus <- GLoutput[GLoutput$value == 0 & !is.na(GLoutput$value), 
        ]
    GLstrata$notest <- GLoutput[is.na(GLoutput$value), ]
    return(GLstrata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
