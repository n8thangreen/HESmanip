\name{survDataByGender}
\alias{survDataByGender}
\docType{data}
\title{Survival analysis format stratified by gender}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(survDataByGender)}
\format{
  The format is:
List of 2 Male and Female.

 A data frame with 100 observations on the following 29 variables.
  \describe{
    \item{\code{hes_disdestdescription}}{a factor with levels \code{} \code{Court} \code{Death, including still birth} \code{High Security psych. Hospital, Scotland} \code{Loc. Auth. foster care, but not in Part 3 res. accommodation where care is provided} \code{Loc. Auth. Part 3 res. accommodation where care is provided} \code{NHS other hospital provider - high security psychiatric accommodation (from 1 April 2004)} \code{NHS other hospital provider - medium secure unit (from 1 April 2004)} \code{NHS other hospital provider - ward for general patients or the younger phys. disabled (from 1 April 2004)} \code{NHS other hospital provider - ward for maternity patients or neonates (from 1 April 2004)} \code{NHS other hospital provider - ward for patients who are mentally ill or have learning disabilities (from 1 April 2004)} \code{NHS run care home (from 1 April 2004)} \code{Non-NHS (other than Loc. Auth.) run care home} \code{Non-NHS (other than Loc. Auth.) run hospice} \code{Non-NHS run hospital} \code{Non-NHS run hospital - medium secure unit} \code{Not applicable} \code{Penal establishment - police station (from 1999-2000)} \code{Repatriation from a High Security psych. Hospital (HSPH)} \code{Temporary place of residence when usually resident elsewhere, for example, hotels and res. education} \code{The usual place of residence, including no fixed abode} \code{Unknown}}
    \item{\code{hes_digitprovidercode}}{a factor with levels \code{5C400} \code{5CN00} \code{5CQ00} \code{5EF} \code{5F100} \code{5FE00} \code{5JE99} \code{5K3} \code{5KG00} \code{5L100} \code{5M121} \code{5M300} \code{5M700} \code{5MD00} \code{5MV00} \code{5N6AH} \code{5N6AJ} \code{5N6AK} \code{5N6AL} \code{5N6AM} \code{5N6AN} \code{5N6NR} \code{5N6RH} \code{5NE} \code{5NV00} \code{5NW00} \code{5P500} \code{5P537} \code{5P6} \code{5P600} \code{5PAAC} \code{5PACR} \code{5PACT} \code{5PACW} \code{5PADE} \code{5PE} \code{5PJ00} \code{5PK} \code{5PL} \code{5PM} \code{5PM00} \code{5PP} \code{5PP00} \code{5PQ00} \code{5QC00} \code{5QD00} \code{5QH00} \code{5QK} \code{5QK00} \code{5QM} \code{5QM00} \code{5QP00} \code{5QQ00} \code{5QT} \code{5QT00} \code{NT301} \code{NT304} \code{NT308} \code{NT314} \code{NT315} \code{NT318} \code{NT319} \code{NT320} \code{NT321} \code{NT322} \code{NT324} \code{NT325} \code{NT332} \code{NT333} \code{NT337} \code{NT339} \code{NT340} \code{NT347} \code{NT348} \code{NT350} \code{NT351} \code{NT403} \code{NT420} \code{NT501} \code{NT502} \code{NT503} \code{NT506} \code{NT507} \code{NT602} \code{NT603} \code{NT604} \code{NT714} \code{NTA01} \code{NTA03} \code{NTA04} \code{NTC01} \code{NTD01} \code{NTD02} \code{NTE03} \code{NTH11} \code{NTY57} \code{NTYE7} \code{NTYH9} \code{NVC01} \code{NVC04} \code{NVC05} \code{NVC06} \code{NVC07} \code{NVC08} \code{NVC10} \code{NVC12} \code{NVC13} \code{NVC14} \code{NVC15} \code{NVC16} \code{NVC17} \code{NVC18} \code{NVC19} \code{NVC20} \code{NVC21} \code{NVC22} \code{NVC23} \code{NVC24} \code{NVC25} \code{NVC27} \code{NVC28} \code{NVC29} \code{NVC31} \code{NVC33} \code{NVC34} \code{NVC35} \code{NVC40} \code{NVG01} \code{RA201} \code{RA300} \code{RA400} \code{RA7} \code{RA700} \code{RA900} \code{RAE00} \code{RAJ01} \code{RAJ12} \code{RAL00} \code{RAL01} \code{RAL16} \code{RALMV} \code{RAN00} \code{RAPNM} \code{RAS00} \code{RAT} \code{RAX01} \code{RBA00} \code{RBB00} \code{RBD01} \code{RBF00} \code{RBK00} \code{RBL} \code{RBL00} \code{RBN01} \code{RBQ00} \code{RBS25} \code{RBT00} \code{RBV01} \code{RBZ00} \code{RC100} \code{RC368} \code{RC900} \code{RCB00} \code{RCC00} \code{RCD00} \code{RCF00} \code{RCUEF} \code{RCX00} \code{RD100} \code{RD300} \code{RD7} \code{RD800} \code{RDDH0} \code{RDDH1} \code{RDDH8} \code{RDEE2} \code{RDEE4} \code{RDEEB} \code{RDEEE} \code{RDENF} \code{RDEP1} \code{RDEP4} \code{RDU01} \code{RDY00} \code{RDZ} \code{RE900} \code{REF00} \code{REM00} \code{REN20} \code{REP01} \code{REP21} \code{RET20} \code{RF4DG} \code{RF4PH} \code{RF4QH} \code{RFF00} \code{RFR00} \code{RFRDR} \code{RFRPA} \code{RFSDA} \code{RFW00} \code{RFW01} \code{RG200} \code{RG3} \code{RG300} \code{RGCKH} \code{RGD00} \code{RGM21} \code{RGN00} \code{RGP} \code{RGP00} \code{RGQ00} \code{RGR00} \code{RGT} \code{RGT01} \code{RGZ} \code{RGZ01} \code{RH500} \code{RH800} \code{RHA} \code{RHM00} \code{RHQCC} \code{RHQHH} \code{RHQNG} \code{RHQWP} \code{RHU00} \code{RHW01} \code{RJ100} \code{RJ224} \code{RJ600} \code{RJ700} \code{RJ800} \code{RJC00} \code{RJD01} \code{RJD13} \code{RJDTC} \code{RJE01} \code{RJF02} \code{RJL00} \code{RJN00} \code{RJR00} \code{RJX00} \code{RJZ01} \code{RJZ03} \code{RJZ33} \code{RK5BC} \code{RK5HP} \code{RK900} \code{RKB00} \code{RKEQ4} \code{RKL} \code{RL100} \code{RL400} \code{RLNDU} \code{RLNGL} \code{RLNGM} \code{RLNGT} \code{RLNQE} \code{RLQ} \code{RLT00} \code{RLU00} \code{RLY00} \code{RM100} \code{RM200} \code{RM300} \code{RM400} \code{RMC00} \code{RMP01} \code{RMY00} \code{RN100} \code{RN300} \code{RN506} \code{RN525} \code{RN5T1} \code{RN700} \code{RNA00} \code{RNH} \code{RNH00} \code{RNJ00} \code{RNL00} \code{RNN} \code{RNQ00} \code{RNS00} \code{RNU} \code{RNU00} \code{RNZ00} \code{RP1} \code{RP400} \code{RP5BA} \code{RP5DR} \code{RP5MM} \code{RP5RH} \code{RP5TR} \code{RP600} \code{RP700} \code{RP7BJ} \code{RPA00} \code{RPAMM} \code{RPAMP} \code{RPC00} \code{RPC04} \code{RPC12} \code{RPG00} \code{RPL00} \code{RPR00} \code{RPY00} \code{RQ300} \code{RQ330} \code{RQ6} \code{RQ601} \code{RQ614} \code{RQ617} \code{RQ800} \code{RQM00} \code{RQQ} \code{RQWG0} \code{RQWG1} \code{RQWG2} \code{RQWG4} \code{RQX} \code{RQX00} \code{RQY00} \code{RR100} \code{RR101} \code{RR105} \code{RR109} \code{RR700} \code{RR800} \code{RRD02} \code{RRE00} \code{RRF00} \code{RRJ00} \code{RRK02} \code{RRK03} \code{RRP00} \code{RRV03} \code{RRV11} \code{RRV13} \code{RRV30} \code{RRV92} \code{RRVE1} \code{RRVNC} \code{RRVNQ} \code{RT100} \code{RT200} \code{RT300} \code{RT502} \code{RT600} \code{RTD00} \code{RTE00} \code{RTF00} \code{RTG} \code{RTH00} \code{RTK00} \code{RTP00} \code{RTQ00} \code{RTR00} \code{RTR45} \code{RTV00} \code{RTX00} \code{RV3} \code{RV504} \code{RV700} \code{RV800} \code{RV900} \code{RVJ00} \code{RVL00} \code{RVL01} \code{RVL07} \code{RVLC7} \code{RVN00} \code{RVR00} \code{RVV01} \code{RVW00} \code{RVY01} \code{RVY02} \code{RW100} \code{RW102} \code{RW103} \code{RW300} \code{RW401} \code{RW500} \code{RW600} \code{RWA00} \code{RWDAA} \code{RWDDA} \code{RWDLA} \code{RWDLP} \code{RWE00} \code{RWF00} \code{RWG00} \code{RWH} \code{RWH04} \code{RWJ00} \code{RWK00} \code{RWN00} \code{RWP00} \code{RWQ} \code{RWR} \code{RWV00} \code{RWW00} \code{RWX00} \code{RWY00} \code{RX1CC} \code{RX1RA} \code{RX2} \code{RX3} \code{RX301} \code{RX302} \code{RX401} \code{RXA00} \code{RXC00} \code{RXE00} \code{RXF00} \code{RXF03} \code{RXF04} \code{RXF05} \code{RXF10} \code{RXG} \code{RXH00} \code{RXK01} \code{RXK02} \code{RXKTC} \code{RXL00} \code{RXM00} \code{RXN00} \code{RXP00} \code{RXPWE} \code{RXQ00} \code{RXR01} \code{RXR02} \code{RXT00} \code{RXV00} \code{RXWAS} \code{RXWAT} \code{RXWMB} \code{RXWMJ} \code{RXWML} \code{RXWMS} \code{RXWMT} \code{RXX10} \code{RXX21} \code{RXXED} \code{RXY00} \code{RXY04} \code{RYG} \code{RYG00} \code{RYJ} \code{RYK} \code{TADLM} \code{TAE05} \code{TAF00} \code{TAH00} \code{TAJ11} \code{TAJ20} \code{RFFAA} \code{RQ600} \code{RVV09} \code{RVV99} \code{RWF03} \code{RYQ10} \code{RYQ30} \code{RYR00}}
    \item{\code{hes_neocaredescription}}{a factor with levels \code{} \code{Level 1 intensive care} \code{Level 2 intensive care} \code{Normal care} \code{Not applicable} \code{Not known} \code{Special care}}
    \item{\code{hes_admimethdescription}}{a factor with levels \code{Elective - booked} \code{Elective - from waiting list} \code{Elective - planned} \code{Emergency - other means, including patients who arrive via A&E department of another HC provider} \code{Emergency - via A&E services, including casualty department of provider} \code{Emergency - via Bed Bureau, including Central Bureau} \code{Emergency - via consultant out-patient clinic} \code{Emergency - via General Practitioner (GP)} \code{Maternity - where baby was delivered after mother's admission} \code{Maternity - where baby was delivered before mother's admission} \code{Not known} \code{Other - babies born in HC provider} \code{Other - babies born outside HC provider, except when born at home as intended} \code{Transfer of any admitted patient from another hospital provider}}
    \item{\code{hes_admisorcdescription}}{a factor with levels \code{Babies born in or on the way to hospital} \code{Loc. Auth. foster care, but not in Part 3 res. accommodation (where care is provided)} \code{Loc. Auth. Part 3 res. accommodation (where care is provided)} \code{NHS other  hospital provider - ward for general patients or the younger phys. disabled or A&E depatment (from 1 April 2004)} \code{NHS other hospital provider - high security psychiatric accommodation (from 1 April 2004)} \code{NHS other hospital provider - ward for maternity patients or neonatal ward (from 1 April 2004)} \code{NHS other hospital provider - ward for patients who are mentally ill or have learning disabilities (from 1 April 2004)} \code{NHS run care home (from 1 April 2004)} \code{Non-NHS (other than Loc. Auth.) run care home} \code{Non-NHS (other than Loc. Auth.) run hospice} \code{Non-NHS run hospital} \code{Not known} \code{Penal establishment, court or police station (Court only from 1999&2000)} \code{Temporary place of residence when usually resident elsewhere} \code{The usual place of residence, including no fixed abode}}
    \item{\code{surgical}}{a logical vector}
    \item{\code{cong}}{a logical vector}
    \item{\code{head}}{a logical vector}
    \item{\code{infection}}{a logical vector}
    \item{\code{cath}}{a logical vector}
    \item{\code{prem}}{a logical vector}
    \item{\code{mental}}{a logical vector}
    \item{\code{other}}{a logical vector}
    \item{\code{Tai}}{a logical vector}
    \item{\code{unknown}}{a logical vector}
    \item{\code{cancer}}{a logical vector}
    \item{\code{highRisk}}{a logical vector}
    \item{\code{res}}{a logical vector}
    \item{\code{infstatus}}{a numeric vector}
    \item{\code{missingType}}{a numeric vector}
    \item{\code{missingTime}}{a logical vector}
    \item{\code{agegr}}{a factor with levels \code{(100,1e+04]} \code{(0,2]} \code{(2,5]} \code{(5,10]} \code{(10,18]} \code{(18,40]} \code{(40,60]} \code{(60,100]}}
    \item{\code{spectime}}{a numeric vector}
    \item{\code{gender}}{a factor with levels \code{Female} \code{Male} \code{Not known} \code{Not specified}}
    \item{\code{age}}{a numeric vector}
    \item{\code{event}}{a logical vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{time2death}}{a numeric vector}
    \item{\code{time2disch}}{a numeric vector}
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
data(survDataByGender)
## maybe str(survDataByGender) ; plot(survDataByGender) ...
}
\keyword{datasets}
