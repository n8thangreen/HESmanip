\name{infdata_sample}
\alias{infdata_sample}
\docType{data}
\title{Exposure/infected sample of patients in linked HES-LabBase2 format}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(infdata_sample)}
\format{
  A data frame with 20 observations on the following 100 variables.
  \describe{
    \item{\code{lab_opieid}}{a factor with levels \code{09ISW034654/01} \code{09LIV003872/01} \code{09LIV006355/01} \code{09LIV009557/01} \code{09LIV016856/01} \code{09LIV020231/01} \code{09SLY017471/01} \code{09ZWT025960/01} \code{09ZXC023311/01} \code{09ZXK039504/01} \code{09ZXK039700/01} \code{09ZXK040794/01} \code{10BOL007094/01} \code{10FRE057722/01} \code{10LIV045064/01} \code{10MMM021350/01} \code{10NOT137976/01} \code{10ZWG152209/01} \code{10ZWG165532/01}}
    \item{\code{hes_ID}}{a numeric vector}
    \item{\code{hes_hesid_count_duplicates}}{a numeric vector}
    \item{\code{lab_organismname}}{a factor with levels \code{COLIFORM} \code{ENTEROBACTER CLOACAE} \code{ENTEROCOCCUS FAECALIS} \code{ENTEROCOCCUS SP} \code{KLEBSIELLA PNEUMONIAE} \code{MORAXELLA CATARRHALIS} \code{PANTOEA SP} \code{PSEUDOMONAS AERUGINOSA} \code{STAPHYLOCOCCUS COAGULASE NEGATIVE}}
    \item{\code{lab_group}}{a factor with levels \code{COAGULASE NEGATIVE STAPHYLOCOCCUS} \code{ENTEROBACTER} \code{ENTEROCOCCUS} \code{KLEBSIELLA} \code{other} \code{P. AERUGINOSA}}
    \item{\code{lab_GRAMEXT}}{a factor with levels \code{N} \code{P}}
    \item{\code{X}}{a numeric vector}
    \item{\code{lab_Specimendate}}{a factor with levels \code{2009-01-09} \code{2009-06-04} \code{2009-08-11} \code{2009-08-25} \code{2009-08-28} \code{2009-09-01} \code{2009-09-15} \code{2009-09-27} \code{2009-10-23} \code{2009-10-29} \code{2009-11-11} \code{2009-11-27} \code{2009-12-08} \code{2010-01-07} \code{2010-01-30} \code{2010-02-08} \code{2010-02-25} \code{2010-03-04}}
    \item{\code{hes_admdte}}{a factor with levels \code{2008-10-02} \code{2009-05-24} \code{2009-07-20} \code{2009-08-18} \code{2009-08-22} \code{2009-09-07} \code{2009-09-23} \code{2009-10-12} \code{2009-11-01} \code{2009-12-05} \code{2010-01-04} \code{2010-01-05} \code{2010-01-21} \code{2010-03-02}}
    \item{\code{hes_disdte}}{a factor with levels \code{2009-03-11} \code{2009-06-16} \code{2009-09-21} \code{2009-09-23} \code{2009-10-13} \code{2009-11-10} \code{2009-11-21} \code{2010-01-10} \code{2010-01-18} \code{2010-01-30} \code{2010-02-16} \code{2010-03-11} \code{2010-03-17}}
    \item{\code{adm_diff_specdate}}{a numeric vector}
    \item{\code{dis_diff_specdate}}{a numeric vector}
    \item{\code{adm_diff_dis}}{a numeric vector}
    \item{\code{lab_AMPAMOX}}{a numeric vector}
    \item{\code{lab_CTX}}{a numeric vector}
    \item{\code{lab_CTZ}}{a numeric vector}
    \item{\code{lab_CTRI}}{a numeric vector}
    \item{\code{lab_GENT}}{a numeric vector}
    \item{\code{lab_IMI}}{a numeric vector}
    \item{\code{lab_MERO}}{a numeric vector}
    \item{\code{lab_MET}}{a numeric vector}
    \item{\code{lab_PIPTAZ}}{a numeric vector}
    \item{\code{lab_labname}}{a factor with levels \code{CHELSEA AND WESTMINSTER HOSPITAL (LONDON)} \code{HILLINGDON HOSPITAL (MIDDLESEX)} \code{HPA NORTH WEST  MANCHESTER LABORATORY} \code{IPSWICH MICROBIOLOGY LABORATORY} \code{NOTTINGHAM MICROBIOLOGY LABORATORY} \code{SALISBURY MICROBIOLOGY LABORATORY} \code{SOUTHMEAD HOSPITAL (BRISTOL)} \code{ST HELIER HOSPITAL (CARSHALTON)} \code{ST MARY'S HOSPITAL (LONDON)} \code{UNIVERSITY COLLEGE HOSPITAL (LONDON)}}
    \item{\code{lab_currentregionname}}{a factor with levels \code{E MIDS} \code{EAST} \code{LONDON} \code{N WEST} \code{S WEST}}
    \item{\code{lab_CoNSmet}}{a numeric vector}
    \item{\code{lab_AMINOAMOX1}}{a numeric vector}
    \item{\code{lab_AMINOMET1}}{a numeric vector}
    \item{\code{lab_AMINOVAN1}}{a numeric vector}
    \item{\code{lab_CTXI1}}{a numeric vector}
    \item{\code{lab_pen}}{a numeric vector}
    \item{\code{lab_met}}{a numeric vector}
    \item{\code{lab_van2}}{a numeric vector}
    \item{\code{lab_AMINOAMOXAMP1}}{a numeric vector}
    \item{\code{lab_AMINOPIPTAZ1}}{a numeric vector}
    \item{\code{lab_AMINOCTZ1}}{a numeric vector}
    \item{\code{lab_AMINOIMI1}}{a numeric vector}
    \item{\code{lab_AMINOMER1}}{a numeric vector}
    \item{\code{lab_CTXI2}}{a numeric vector}
    \item{\code{hes_sexdescription}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{hes_ageatstartofepisode}}{a numeric vector}
    \item{\code{hes_neocaredescription}}{a factor with levels \code{Not applicable} \code{Not known} \code{Special care}}
    \item{\code{hes_admimethdescription}}{a factor with levels \code{Elective - booked} \code{Elective - planned} \code{Emergency - other means, including patients who arrive via A&E department of another HC provider} \code{Maternity - where baby was delivered after mother's admission} \code{Other - babies born in HC provider} \code{Transfer of any admitted patient from another hospital provider}}
    \item{\code{hes_admisorcdescription}}{a factor with levels \code{Babies born in or on the way to hospital} \code{NHS other  hospital provider - ward for general patients or the younger phys. disabled or A&E depatment (from 1 April 2004)} \code{NHS other hospital provider - ward for maternity patients or neonatal ward (from 1 April 2004)} \code{The usual place of residence, including no fixed abode}}
    \item{\code{hes_dismethdescription}}{a factor with levels \code{Died} \code{Discharged on clinical advice or with clinical consent}}
    \item{\code{hes_dischargeddeadoralive}}{a factor with levels \code{Death} \code{Live}}
    \item{\code{X.1}}{a numeric vector}
    \item{\code{hes_disdestdescription}}{a factor with levels \code{Death, including still birth} \code{NHS other hospital provider - ward for general patients or the younger phys. disabled (from 1 April 2004)} \code{NHS other hospital provider - ward for maternity patients or neonates (from 1 April 2004)} \code{NHS run care home (from 1 April 2004)} \code{The usual place of residence, including no fixed abode}}
    \item{\code{hes_diag1description}}{a factor with levels \code{Acute lymphoblastic leukaemia} \code{Acute myeloid leukaemia} \code{Adenoviral pneumonia} \code{Congenital absence, atresia and stenosis of ileum} \code{Crohn's disease of small intestine} \code{Extremely low birth weight} \code{Labour and delivery complicated by fetal heart rate anomaly} \code{Malignant neoplasm, other sarcomas of liver} \code{Other and unspecified intestinal obstruction} \code{Other low birth weight} \code{Other specified congenital malformations of skin} \code{Postoperative intestinal obstruction}}
    \item{\code{hes_diag1_code_primary}}{a factor with levels \code{C224-} \code{C910-} \code{C920} \code{J120-} \code{K500-} \code{K566} \code{K913-} \code{O680} \code{P070} \code{P071} \code{Q412} \code{Q828}}
    \item{\code{hes_diag2description}}{a factor with levels \code{} \code{Abnormal weight loss} \code{Agranulocytosis} \code{Anaemia, unspecified} \code{Colostomy and enterostomy malfunction} \code{Down's syndrome, unspecified} \code{Hepatic fibrosis} \code{Hyperkinetic disorder, unspecified} \code{Intestinal adhesions [bands] with obstruction} \code{Necrotizing enterocolitis of fetus and newborn} \code{Pleural effusion, not elsewhere classified} \code{Respiratory syncytial virus pneumonia} \code{Single live birth} \code{Staph aureus as cause of dis classified to other chapters} \code{Twin, born in hospital}}
    \item{\code{hes_diag2_code}}{a factor with levels \code{} \code{B956} \code{D649-} \code{D70X-} \code{F909-} \code{J121-} \code{J90X} \code{K565-} \code{K740} \code{K914-} \code{P77X} \code{Q909-} \code{R634} \code{Z370} \code{Z383}}
    \item{\code{hes_diag3description}}{a factor with levels \code{} \code{Acute bronchiolitis due to respiratory syncytial virus} \code{Agranulocytosis} \code{Anaemia, unspecified} \code{Bronchopulmonary dysplasia origin in the perinatal period} \code{Childhood autism} \code{Essential (primary) hypertension} \code{Flatulence and related conditions} \code{Glossodynia} \code{Hyperkinetic disorder, unspecified} \code{Infection following a procedure, not elsewhere classified} \code{Necrotizing enterocolitis of fetus and newborn} \code{Other generalized epilepsy and epileptic syndromes} \code{Pyrexia of unknown origin following delivery} \code{Sepsis of newborn due to Escherichia coli}}
    \item{\code{hes_diag3_code}}{a factor with levels \code{} \code{D649} \code{D70X-} \code{F840-} \code{F909-} \code{G404-} \code{I10X} \code{J210} \code{K146-} \code{O864} \code{P271} \code{P364} \code{P77X-} \code{R14X} \code{T814-}}
    \item{\code{hes_diag4description}}{a factor with levels \code{} \code{Childhood autism} \code{Developmental disorder of scholastic skills, unspecified} \code{Down's syndrome, unspecified} \code{Fever, unspecified} \code{Fluid overload} \code{Haemorrhage into co-twin} \code{Other bacterial meningitis} \code{Other preterm infants} \code{Other spec bact agents as cause of dis class oth chaps} \code{Other spec dis and conditions comp preg childbirth puerp} \code{Other specified septicaemia} \code{Primary pulmonary hypertension} \code{Pulmonary collapse} \code{Thrombocytopenia, unspecified}}
    \item{\code{hes_diag4_code}}{a factor with levels \code{} \code{A418} \code{B968-} \code{D696} \code{E877-} \code{F819} \code{F840-} \code{G008} \code{I270-} \code{J981-} \code{O998} \code{P073-} \code{P503} \code{Q909-} \code{R509}}
    \item{\code{hes_diag5description}}{a factor with levels \code{} \code{Chest pain, unspecified} \code{Decubitus ulcer} \code{Dyspnoea} \code{Escherichia coli as cause of dis classified to other chaps} \code{Fet and newborn affected by placental transfusion syndr} \code{Fever, unspecified} \code{Gastro-oesophageal reflux disease without oesophagitis} \code{Hypotension, unspecified} \code{Noninfective gastroenteritis and colitis, unspecified} \code{Ocular pain} \code{Personal hist of cong malf'ns deform'ns/chromosomal abnms} \code{Thrombocytopenia, unspecified}}
    \item{\code{hes_diag5_code}}{a factor with levels \code{} \code{B962} \code{D696} \code{H571-} \code{I959-} \code{K219-} \code{K529-} \code{L89X} \code{P023} \code{R060} \code{R074-} \code{R509-} \code{Z877}}
    \item{\code{hes_diag6description}}{a factor with levels \code{} \code{Acquired absence of other parts of digestive tract} \code{Agranulocytosis} \code{Chemotherapy session for neoplasm} \code{Haematemesis} \code{Hypo-osmolality and hyponatraemia} \code{Low back pain} \code{Neonatal jaundice associated with preterm delivery} \code{Other disorders of bilirubin metabolism} \code{Other feeding problems of newborn} \code{Other lack of expected normal physiological development}}
    \item{\code{hes_diag6_code}}{a factor with levels \code{} \code{D70X} \code{E806} \code{E871-} \code{K920-} \code{M5458} \code{P590} \code{P928} \code{R628-} \code{Z511-} \code{Z904-}}
    \item{\code{hes_diag7description}}{a factor with levels \code{} \code{Epistaxis} \code{Melaena} \code{Other antineopastic drugs} \code{Other lack of expected normal physiological development} \code{Other viral infections of unspecified site} \code{Personal history of allergy to analgesic agent} \code{Personal history of diseases of the digestive system}}
    \item{\code{hes_diag7_code}}{a factor with levels \code{} \code{B348-} \code{K921-} \code{R040} \code{R628-} \code{Y433} \code{Z871-} \code{Z886-}}
    \item{\code{hes_diag8description}}{a factor with levels \code{} \code{Anaemia, unspecified} \code{Cough} \code{Infect inflamm reac due oth card vasc devs implant and graft} \code{Infection following a procedure, not elsewhere classified} \code{Personal history of dis of the nervous sys and sense organs}}
    \item{\code{hes_diag8_code}}{a factor with levels \code{} \code{D649} \code{R05X-} \code{T814-} \code{T827} \code{Z866-}}
    \item{\code{hes_diag9description}}{a factor with levels \code{} \code{Hepatomegaly with splenomegaly, not elsewhere classified} \code{Klebsiella pneumoniae as cause dis class other chaps} \code{Oedema, unspecified} \code{Other and unspecified abdominal pain}}
    \item{\code{hes_diag9_code}}{a factor with levels \code{} \code{B961-} \code{R104-} \code{R162} \code{R609}}
    \item{\code{hes_diag10description}}{a factor with levels \code{} \code{Acute nasopharyngitis [common cold]} \code{Megacolon, not elsewhere classified} \code{Nausea and vomiting} \code{Septicaemia due to other Gram-negative organisms}}
    \item{\code{hes_diag10_code}}{a factor with levels \code{} \code{A415-} \code{J00X} \code{K593} \code{R11X-}}
    \item{\code{hes_diag11description}}{a factor with levels \code{} \code{Fever, unspecified} \code{Hypokalaemia}}
    \item{\code{hes_diag11_code}}{a factor with levels \code{} \code{E876} \code{R509-}}
    \item{\code{hes_diag12description}}{a factor with levels \code{} \code{Pericardial effusion (noninflammatory)} \code{Tachycardia, unspecified}}
    \item{\code{hes_diag12_code}}{a factor with levels \code{} \code{I313-} \code{R000}}
    \item{\code{hes_diag13description}}{a factor with levels \code{} \code{Anxiety disorder, unspecified} \code{Carrier of other specified bacterial diseases}}
    \item{\code{hes_diag13_code}}{a factor with levels \code{} \code{F419} \code{Z223-}}
    \item{\code{hes_diag14description}}{a logical vector}
    \item{\code{hes_diag14_code}}{a logical vector}
    \item{\code{hes_diag15description}}{a logical vector}
    \item{\code{hes_diag15_code}}{a logical vector}
    \item{\code{hes_diag16description}}{a logical vector}
    \item{\code{hes_diag16_code}}{a logical vector}
    \item{\code{hes_diag17description}}{a logical vector}
    \item{\code{hes_diag17_code}}{a logical vector}
    \item{\code{hes_diag18description}}{a logical vector}
    \item{\code{hes_diag18_code}}{a logical vector}
    \item{\code{hes_diag19description}}{a logical vector}
    \item{\code{hes_diag19_code}}{a logical vector}
    \item{\code{hes_diag20description}}{a logical vector}
    \item{\code{hes_diag20_code}}{a logical vector}
    \item{\code{hes_oper1description}}{a factor with levels \code{-} \code{All normal delivery} \code{Delivery of subsequent element of cycle of chemotherapy for neoplasm} \code{Ileectomy and anastomosis of ileum to ileum} \code{Injection of therapeutic substance into cerebrospinal fluid} \code{Intravenous nutrition Band 1} \code{Non-invasive ventilation NEC} \code{Resection of ileo-colic anastomosis} \code{Unspecified diagnostic spinal puncture}}
    \item{\code{hes_op1_maincode}}{a factor with levels \code{-} \code{A542} \code{A559} \code{E852} \code{G693} \code{G734} \code{R249} \code{X724} \code{X904}}
    \item{\code{hes_oper2description}}{a factor with levels \code{} \code{Diagnostic extraction of bone marrow NEC} \code{Episiotomy to facilitate delivery} \code{Freeing of adhesions of peritoneum} \code{Intravenous immunotherapy} \code{Intravenous nutrition Band 1} \code{Lumbar spinal cord} \code{Other specified procurement of drugs for chemotherapy for neoplasm in Bands 1-5}}
    \item{\code{hes_op2_code}}{a factor with levels \code{} \code{R271} \code{T413} \code{W365} \code{X353} \code{X708} \code{X904} \code{Z063}}
    \item{\code{hes_oper3description}}{a factor with levels \code{} \code{Intravenous nutrition Band 1} \code{Other specified procurement of drugs for chemotherapy for neoplasm in Bands 1-5} \code{Transthoracic echocardiography} \code{Unspecified creation of artificial opening into ileum} \code{Unspecified phototherapy to skin} \code{Wing of ilium}}
    \item{\code{hes_op3_code}}{a factor with levels \code{} \code{G749} \code{S129} \code{U201} \code{X708} \code{X904} \code{Z753}}
    \item{\code{hes_oper4description}}{a factor with levels \code{} \code{Cytokine inhibitor drugs Band 1} \code{Intravenous blood transfusion of packed cells} \code{Neutropenia drugs Band 1} \code{Open insertion of central venous catheter} \code{Right sided operation}}
    \item{\code{hes_op4_code}}{a factor with levels \code{} \code{L911} \code{X332} \code{X903} \code{X921} \code{Z942}}
    \item{\code{hes_oper5description}}{a factor with levels \code{} \code{Antifungal drugs Band 2} \code{Attention to central venous catheter NEC} \code{Procurement of drugs for chemotherapy for neoplasm for regimens in Band 1} \code{Specified vein NEC}}
    \item{\code{hes_op5_code}}{a factor with levels \code{} \code{L913} \code{X701} \code{X862} \code{Z398}}
    \item{\code{hes_procode5description}}{a factor with levels \code{CENTRAL MANCHESTER UNIVERSITY HOSPITALS NHS FOUNDATION TRUST} \code{CHELSEA AND WESTMINSTER HOSPITAL NHS FOUNDATION TRUST} \code{EPSOM AND ST HELIER UNIVERSITY HOSPITALS NHS TRUST} \code{IMPERIAL COLLEGE HEALTHCARE NHS TRUST} \code{IPSWICH HOSPITAL NHS TRUST} \code{NORTH BRISTOL NHS TRUST} \code{NOTTINGHAM UNIVERSITY HOSPITALS NHS TRUST - QUEEN'S MEDICAL CENTRE CAMPUS} \code{SALISBURY HEALTH CARE NHS TRUST} \code{THE HILLINGDON HOSPITAL NHS TRUST} \code{UNIVERSITY COLLEGE HOSPITAL}}
    \item{\code{hes_digitprovidercode}}{a factor with levels \code{RAS00} \code{RGQ00} \code{RNZ00} \code{RQM00} \code{RRV03} \code{RVJ00} \code{RVR00} \code{RW300} \code{RX1RA} \code{RYJ}}
    \item{\code{hes_admincat_description}}{a factor with levels \code{NHS patient}}
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
data(infdata_sample)
## maybe str(infdata_sample) ; plot(infdata_sample) ...
}
\keyword{datasets}
