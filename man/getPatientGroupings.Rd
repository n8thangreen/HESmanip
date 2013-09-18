\name{getPatientGroupings}
\alias{getPatientGroupings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Group patients by risk factors}
\description{
Use a combination of fields identified through partial word-matching and fields explicitly identified by subject matter experts and look-up table sources.
}
\usage{
getPatientGroupings(inf.data, mix.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf.data}{
data.frame of infected patient data.
}
  \item{mix.data}{
data.frame of non-infected patient data.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
data.frame of logical indicators whether a patient belongs to a risk group.

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
function (inf.data, mix.data) 
{
    SMEcodes <- read.csv(".\\Reference_tables\\explicitSMEcodes.csv")
    codeColnames <- c("Code", "Group", "Description", "Scheme")
    SSI_OPCS_codes <- read.csv(".\\Reference_tables\\SSI_OPCS_codes.csv")
    colnames(SSI_OPCS_codes) <- codeColnames
    SSI_OPCS_codes$Code <- clean.ref(SSI_OPCS_codes$Code)
    opcode_cols.inf <- grep("hes_op([0-9]+)_(main)?code", names(inf.data))
    opcode_cols.mix <- grep("operativeprocedure", names(mix.data))
    inf.data <- clean.data(inf.data, opcode_cols.inf)
    mix.data <- clean.data(mix.data, opcode_cols.mix)
    surgFlags.inf <- facFlags(SSI_OPCS_codes, "OPCS", inf.data[, 
        opcode_cols.inf])
    surgFlags.mix <- facFlags(SSI_OPCS_codes, "OPCS", mix.data[, 
        opcode_cols.mix])
    SMEcodesFlags.inf <- facFlags(SMEcodes, "OPCS", inf.data[, 
        opcode_cols.inf])
    SMEcodesFlags.mix <- facFlags(SMEcodes, "OPCS", mix.data[, 
        opcode_cols.mix])
    diagcode_cols.mix <- seq(13, 51, by = 2)
    diagcode_cols.inf <- grep("hes_diag([0-9]+)_code", names(inf.data))
    inf.data <- clean.data(inf.data, diagcode_cols.inf)
    mix.data <- clean.data(mix.data, diagcode_cols.mix)
    SMEcodesFlags.inf <- SMEcodesFlags.inf | facFlags(SMEcodes, 
        "ICD10", inf.data[, diagcode_cols.inf])
    SMEcodesFlags.mix <- SMEcodesFlags.mix | facFlags(SMEcodes, 
        "ICD10", mix.data[, diagcode_cols.mix])
    ICDGroup.match <- read.table(".\\Reference_tables\\ICDGroup_match.txt")
    ICDgroupFlag.mix <- facFlags(ICDGroup.match, "ICD10", mix.data[, 
        diagcode_cols.mix])
    ICDgroupFlag.inf <- facFlags(ICDGroup.match, "ICD10", inf.data[, 
        diagcode_cols.inf])
    cathCodesgroupFlags.inf <- SMEcodesFlags.inf[, "invasive"]
    cathCodesgroupFlags.mix <- SMEcodesFlags.mix[, "invasive"]
    premCodesgroupFlags.inf <- SMEcodesFlags.inf[, "low birth"]
    premCodesgroupFlags.mix <- SMEcodesFlags.mix[, "low birth"]
    codesIndicators.inf <- data.frame(surg = surgFlags.inf, cong = SMEcodesFlags.inf[, 
        "congenital"], head = SMEcodesFlags.inf[, "head"], infection = SMEcodesFlags.inf[, 
        "infection"], cath = cathCodesgroupFlags.inf, prem = premCodesgroupFlags.inf, 
        mental = SMEcodesFlags.inf[, "mental"], other = SMEcodesFlags.inf[, 
            "other"], Tai = SMEcodesFlags.inf[, "Tai"], unknown = SMEcodesFlags.inf[, 
            "unknown"], cancer = ICDgroupFlag.inf[, "cancer"])
    codesIndicators.mix <- data.frame(surg = surgFlags.mix, cong = SMEcodesFlags.mix[, 
        "congenital"], head = SMEcodesFlags.mix[, "head"], infection = SMEcodesFlags.mix[, 
        "infection"], cath = cathCodesgroupFlags.mix, prem = premCodesgroupFlags.mix, 
        mental = SMEcodesFlags.mix[, "mental"], other = SMEcodesFlags.mix[, 
            "other"], Tai = SMEcodesFlags.mix[, "Tai"], unknown = SMEcodesFlags.mix[, 
            "unknown"], cancer = ICDgroupFlag.mix[, "cancer"])
    codesIndicators.inf <- cbind(codesIndicators.inf, highRisk = rowSums(codesIndicators.inf) > 
        0)
    codesIndicators.mix <- cbind(codesIndicators.mix, highRisk = rowSums(codesIndicators.mix) > 
        0)
    list(inf = codesIndicators.inf, mix = codesIndicators.mix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
