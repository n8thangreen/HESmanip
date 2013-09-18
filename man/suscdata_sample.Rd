\name{suscdata_sample}
\alias{suscdata_sample}
\docType{data}
\title{Unexposed/uninfected sample of patients in HES format}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(suscdata_sample)}
\format{
  A data frame with 20 observations on the following 64 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{sexdescription}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{ageatstartofepisode}}{a numeric vector}
    \item{\code{neocaredescription}}{a factor with levels \code{Level 1 intensive care} \code{Level 2 intensive care} \code{Normal care} \code{Not applicable} \code{Not known} \code{Special care}}
    \item{\code{admimethdescription}}{a factor with levels \code{Elective - from waiting list} \code{Emergency - other means, including patients who arrive via A&E department of another HC provider} \code{Emergency - via A&E services, including casualty department of provider} \code{Emergency - via General Practitioner (GP)} \code{Other - babies born in HC provider} \code{Transfer of any admitted patient from another hospital provider}}
    \item{\code{admisorcdescription}}{a factor with levels \code{Babies born in or on the way to hospital} \code{NHS other  hospital provider - ward for general patients or the younger phys. disabled or A&E depatment (from 1 April 2004)} \code{The usual place of residence, including no fixed abode}}
    \item{\code{dismethdescription}}{a factor with levels \code{Discharged on clinical advice or with clinical consent} \code{Not applicable patient still in hospital}}
    \item{\code{dischargeddeadoralive}}{a factor with levels \code{} \code{Live}}
    \item{\code{disdestdescription}}{a factor with levels \code{NHS other hospital provider - ward for general patients or the younger phys. disabled (from 1 April 2004)} \code{NHS other hospital provider - ward for maternity patients or neonates (from 1 April 2004)} \code{Not applicable} \code{The usual place of residence, including no fixed abode}}
    \item{\code{admissiondate}}{a factor with levels \code{03/07/2008} \code{05/07/2008} \code{07/08/2008} \code{08/03/2009} \code{11/01/2009} \code{11/03/2009} \code{11/07/2008} \code{12/07/2008} \code{13/02/2009} \code{13/08/2008} \code{14/08/2008} \code{17/07/2008} \code{21/03/2009} \code{24/07/2008} \code{26/07/2008} \code{27/07/2008} \code{30/03/2009}}
    \item{\code{dischargedate}}{a factor with levels \code{} \code{05/07/2008} \code{08/08/2008} \code{10/03/2009} \code{11/03/2009} \code{12/01/2009} \code{12/07/2008} \code{13/02/2009} \code{13/08/2008} \code{16/08/2008} \code{17/07/2008} \code{21/03/2009} \code{21/07/2008} \code{24/07/2008} \code{26/07/2008} \code{27/07/2008} \code{30/03/2009} \code{30/07/2008}}
    \item{\code{diag1description}}{a factor with levels \code{Acute upper respiratory infection, unspecified} \code{Bacterial sepsis of newborn, unspecified} \code{Fever, unspecified} \code{Gastro-oesophageal reflux disease without oesophagitis} \code{Health supervision and care other healthy infant and child} \code{Hirschsprung's disease} \code{Neonatal erythema toxicum} \code{Neonatal skin infection} \code{Other congenital malformations of eyelid} \code{Other spec conditions originating in the perinatal period} \code{Other specified congenital malformations of intestine} \code{Viral intestinal infection, unspecified} \code{Wheezing}}
    \item{\code{primarydiagnosis}}{a factor with levels \code{A084} \code{J069} \code{J069-} \code{K219} \code{K219-} \code{P369-} \code{P394-} \code{P831-} \code{P968} \code{P968-} \code{Q103} \code{Q431} \code{Q438} \code{R062-} \code{R509-} \code{Z762-}}
    \item{\code{diag2description}}{a factor with levels \code{} \code{Congenital adrenogenital disorders associated enzyme def} \code{Neonatal cardiac dysrhythmia} \code{Neonatal jaundice, unspecified} \code{Other atopic dermatitis} \code{Rash and other nonspecific skin eruption} \code{Unspec staphy as cause of dis classif to other chapters} \code{Vomiting in newborn}}
    \item{\code{stsecondarysubsidiarydiagnosis}}{a factor with levels \code{} \code{B958-} \code{E250-} \code{L208} \code{P291} \code{P599-} \code{P920} \code{R21X}}
    \item{\code{diag3description}}{a factor with levels \code{} \code{Other neonatal hypoglycaemia} \code{Procedure not carried out for other reasons} \code{Respiratory distress syndrome of newborn}}
    \item{\code{ndsecondarydiagnosis}}{a factor with levels \code{} \code{P220-} \code{P704-} \code{Z538}}
    \item{\code{diag4description}}{a factor with levels \code{} \code{Abnormal weight loss} \code{Other specified respiratory conditions of newborn}}
    \item{\code{rdsecondarydiagnosis}}{a factor with levels \code{} \code{P288-} \code{R634-}}
    \item{\code{diag5description}}{a factor with levels \code{} \code{Neonatal cardiac dysrhythmia}}
    \item{\code{thsecondarydiagnosis}}{a factor with levels \code{} \code{P291-}}
    \item{\code{diag6description}}{a factor with levels \code{} \code{Disturbances of sodium balance of newborn}}
    \item{\code{v23}}{a factor with levels \code{} \code{P742-}}
    \item{\code{diag7description}}{a factor with levels \code{} \code{Disturbances of potassium balance of newborn}}
    \item{\code{v25}}{a factor with levels \code{} \code{P743-}}
    \item{\code{diag8description}}{a factor with levels \code{} \code{Vomiting in newborn}}
    \item{\code{v27}}{a factor with levels \code{} \code{P920-}}
    \item{\code{diag9description}}{a factor with levels \code{} \code{Other and unspecified convulsions}}
    \item{\code{v29}}{a factor with levels \code{} \code{R568-}}
    \item{\code{diag10description}}{a factor with levels \code{} \code{Bacterial sepsis of newborn, unspecified}}
    \item{\code{v31}}{a factor with levels \code{} \code{P369-}}
    \item{\code{diag11description}}{a logical vector}
    \item{\code{v33}}{a logical vector}
    \item{\code{diag12description}}{a logical vector}
    \item{\code{v35}}{a logical vector}
    \item{\code{diag13description}}{a logical vector}
    \item{\code{v37}}{a logical vector}
    \item{\code{diag14description}}{a logical vector}
    \item{\code{v39}}{a logical vector}
    \item{\code{diag15description}}{a logical vector}
    \item{\code{v41}}{a logical vector}
    \item{\code{diag16description}}{a logical vector}
    \item{\code{v43}}{a logical vector}
    \item{\code{diag17description}}{a logical vector}
    \item{\code{v45}}{a logical vector}
    \item{\code{diag18description}}{a logical vector}
    \item{\code{v47}}{a logical vector}
    \item{\code{diag19description}}{a logical vector}
    \item{\code{v49}}{a logical vector}
    \item{\code{diag20description}}{a logical vector}
    \item{\code{v51}}{a logical vector}
    \item{\code{oper1description}}{a factor with levels \code{-} \code{Other specified continuous Infusion of therapeutic substance} \code{Unspecified diagnostic spinal puncture}}
    \item{\code{mainoperativeprocedure}}{a factor with levels \code{-} \code{A559} \code{X298}}
    \item{\code{oper2description}}{a logical vector}
    \item{\code{stsecondaryoperativeprocedure}}{a logical vector}
    \item{\code{oper3description}}{a logical vector}
    \item{\code{ndsecondaryoperativeprocedure}}{a logical vector}
    \item{\code{oper4description}}{a logical vector}
    \item{\code{rdsecondaryoperativeprocedure}}{a logical vector}
    \item{\code{oper5description}}{a logical vector}
    \item{\code{thsecondaryoperativeprocedure}}{a logical vector}
    \item{\code{procode5description}}{a factor with levels \code{DONCASTER ROYAL INFIRMARY} \code{EAST SUSSEX HOSPITALS NHS TRUST} \code{HEATHERWOOD AND WEXHAM PARK HOSPITALS NHS FOUNDATION TRUST} \code{OXFORD RADCLIFFE HOSPITALS NHS TRUST} \code{PORTSMOUTH HOSPITALS NHS TRUST} \code{SOUTHAMPTON UNIVERSITY HOSPITALS NHS TRUST} \code{ST GEORGE'S HEALTHCARE NHS TRUST} \code{WORTHING AND SOUTHLANDS HOSPITALS NHS TRUST}}
    \item{\code{digitprovidercode}}{a factor with levels \code{RD7} \code{RHM00} \code{RHU00} \code{RJ700} \code{RP5DR} \code{RPL00} \code{RTH00} \code{RXC00}}
    \item{\code{administrativecategory}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(suscdata_sample)
## maybe str(suscdata_sample) ; plot(suscdata_sample) ...
}
\keyword{datasets}
