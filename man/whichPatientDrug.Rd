\name{whichPatientDrug}
\alias{whichPatientDrug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Match patients against a look-up table of treatments}
\description{
Determine which patients are susceptible or resistant to the guideline drugs, accounting for age category and other user defined covariates.
}
\usage{
whichPatientDrug(inf.data, drugGuideline.data, patientdrugs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf.data}{
Linked HES-LabBase2 data set of exposed patients.
}
  \item{drugGuideline.data}{
A data.frame of which antibiotics are first-line treatments for which pathogen. 
}
  \item{patientdrugs}{
List of drugs per patient.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns: A list of susceptible or resistant status per patient: 
  
  \item{Per drug}{TRUE entry if guideline antibiotic and resistant/susceptible}
  \item{Aggregated across all drugs}{TRUE entry if guideline antibiotic and resistant/susceptible}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
N. Green
}
\note{
This is an alternative to the \code{whichPatientDrug.join} function which uses joining operations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AgeCatModify}}, \code{\link{whichPatientDrug.joins}}
}
\examples{
data(infdata_sample)
data(antibioticLookup)

drug.names <<- factor(c("AMPAMOX","CTX","CTZ","CTRI","GENT","IMI","MERO",
                        "MET","PIPTAZ","CoNSmet","AMINOAMOX1","AMINOMET1","AMINOVAN1","CTXI1",
                        "pen","met","van2","AMINOAMOXAMP1","AMINOPIPTAZ1","AMINOCTZ1","AMINOIMI1","AMINOMER1","CTXI2"))

##################
# 0: susceptible #
# 1: resistant   #
# NA: not tested #
##################

## extract drug columns from original inf.data dataset
inf.drug.data <- infdata_sample[, paste("lab_",drug.names,sep="")]
inf.drug.data <- data.frame(lapply(inf.drug.data, factor))

## make a list of the drugs administered to each patient
## for each response type

patientdrugs.res <- patientdrugs.sus <- list()
for (i in 1:nrow(inf.drug.data)){ # each patient
  patientdrugs.res[[i]] <- drug.names[!is.na(inf.drug.data[i,]) & inf.drug.data[i,]==1]  # resistant 
  patientdrugs.sus[[i]] <- drug.names[!is.na(inf.drug.data[i,]) & inf.drug.data[i,]==0]  # susceptible
}

whichPatientDrug.res <- whichPatientDrug(infdata_sample, antibioticLookup.dat, patientdrugs.res)
whichPatientDrug.sus <- whichPatientDrug(infdata_sample, antibioticLookup.dat, patientdrugs.sus)





## The function is currently defined as
function (inf.data, drugGuideline.data, patientdrugs) 
{
    npatient <- length(patientdrugs)
    ndrugs <- length(drug.names)
    nguidelines <- nrow(drugGuideline.data)
    FollowGuideline.tot <- rep(FALSE, npatient)
    FollowGuideline.indiv <- matrix(NA, nrow = npatient, ncol = ndrugs)
    colnames(FollowGuideline.indiv) <- drug.names
    patient.age_cat <- AgeCatModify(inf.data)
    for (i in 1:nguidelines) {
        if (drugGuideline.data$drug_code[i] != "" & drugGuideline.data$drug_code[i] \%in\% 
            drug.names) {
            drugMatch <- sapply(patientdrugs, function(x) drugGuideline.data$drug_code[i] \%in\% 
                x) & (drugGuideline.data$ageCat[i] == patient.age_cat | 
                drugGuideline.data$ageCat[i] == "") & (toupper(drugGuideline.data$suspected_organism[i]) == 
                inf.data$lab_organismname | toupper(drugGuideline.data$suspected_organism[i]) == 
                inf.data$lab_group | toupper(drugGuideline.data$suspected_organism[i]) == 
                "")
            drugCol <- as.character(drugGuideline.data$drug_code[i])
            FollowGuideline.indiv[, drugCol] <- FollowGuideline.indiv[, 
                drugCol] | drugMatch
            FollowGuideline.tot <- FollowGuideline.tot | drugMatch
        }
    }
    list(indiv = FollowGuideline.indiv, all = FollowGuideline.tot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
