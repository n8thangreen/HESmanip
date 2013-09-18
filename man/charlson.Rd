\name{charlson}
\alias{charlson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Charlson comorbidity score
}
\description{
Group patient records using ICD10 (or ICD9) codes and then take weighted sum to give a single comorbidity measure.
}
\usage{
charlson(Patient, icd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Patient}{Patient ID could be Spell ID, Trust ID or HES ID for example}
  \item{icd}{vector of ICD codes}
}
\details{
Three digit ICD codes (e.g. I50) can be matched as well as 4 digits (e.g. I500).
}
\value{
A matrix of unique patient ID by individual Charlson group counts and overall score.
}
\references{
V. Sundararajan et al. / Journal of Clinical Epidemiology 57 (2004) 1288–1294
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
function (Patient, icd) 
{
    d <- cbind.data.frame(Patient, icd)
    nobs <- nrow(d)
    d$code <- rep(NA, nobs)
    icd10 <- list()
    ICD3to4 <- function(x) union(unlist(sapply(x, function(y) if (nchar(y) == 
        3) {
        paste(y, 0:9, sep = "")
    })), x)
    icd10$codes$mi <- c("I21", "I22", "I252")
    icd10$codes$chf <- c("I50")
    icd10$codes$pvd <- c("I71", "I790", "I739", "R02", "Z958", 
        "Z959")
    icd10$codes$cvd <- c("I60", "I61", "I62", "I63", "I65", "I66", 
        "G450", "G451", "G452", "G458", "G459", "G46", "I64", 
        "G454", "I670", "I671", "I672", "I674", "I675", "I676", 
        "I677", "I678", "I679", "I681", "I682", "I688", "I69")
    icd10$codes$dem <- c("F00", "F01", "F02", "F051")
    icd10$codes$cpd <- c("J40", "J41", "J42", "J44", "J43", "J45", 
        "J46", "J47", "J67", "J44", "J60", "J61", "J62", "J63", 
        "J66", "J64", "J65")
    icd10$codes$rhm <- c("M32", "M34", "M332", "M053", "M058", 
        "M059", "M060", "M063", "M069", "M050", "M052", "M051", 
        "M353")
    icd10$codes$pep <- c("K25", "K26", "K27", "K28")
    icd10$codes$mld <- c("K702", "K703", "K73", "K717", "K740", 
        "K742", "K746", "K743", "K744", "K745")
    icd10$codes$dnc <- c("E109", "E119", "E139", "E149", "E101", 
        "E111", "E131", "E141", "E105", "E115", "E135", "E145")
    icd10$codes$dwc <- c("E102", "E112", "E132", "E142", "E103", 
        "E113", "E133", "E143", "E104", "E114", "E134", "E144")
    icd10$codes$ple <- c("G81", "G041", "G820", "G821", "G822")
    icd10$codes$ren <- c("N03", "N052", "N053", "N054", "N055", 
        "N056", "N072", "N073", "N074", "N01", "N18", "N19", 
        "N25")
    icd10$codes$can <- c("C0", "C1", "C2", "C3", "C40", "C41", 
        "C43", "C45", "C46", "C47", "C48", "C49", "C5", "C6", 
        "C70", "C71", "C72", "C73", "C74", "C75", "C76", "C80", 
        "C81", "C82", "C83", "C84", "C85", "C883", "C887", "C889", 
        "C900", "C901", "C91", "C92", "C93", "C940", "C941", 
        "C942", "C943", "C9451", "C947", "C95", "C96")
    icd10$codes$met <- c("C77", "C78", "C79", "C80")
    icd10$codes$liv <- c("K729", "K766", "K767", "K721")
    icd10$codes$hiv <- c("B20", "B21", "B22", "B23", "B24")
    icd10$weights <- split(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 
        2, 2, 2, 3, 3, 6), names(icd10$codes))
    icd10$codes <- lapply(icd10$codes, ICD3to4)
    for (i in 1:nobs) {
        patientGroup <- unlist(lapply(icd10$codes, function(x) d[i, 
            "icd"] \%in\% x))
        d$code[i] <- ifelse(sum(patientGroup) == 0, NA, names(icd10$codes)[patientGroup])
    }
    mydata <- cast(d, value = "icd", Patient ~ code, fun.aggregate = length)
    chCols <- mydata[, -c(1, ncol(mydata))] > 0
    comorbidity.n <- rowSums(chCols)
    mydata$charlson <- apply((1 * chCols), 1, function(x) x \%*\% 
        unlist(icd10$weights)[colnames(chCols)])
    mydata$comorbidity.n <- comorbidity.n
    return(mydata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
