charlson <- function (Patient, icd) 
{
  ##
  ##
  ## write.csv(charlson.mix, file=".\\output\\charlson_mix.csv")
  ## charlson.mix <- charlson(Patient=survData.mix$Spell, icd=survData.mix$Prim.Diag.Code)
  
  d <- cbind.data.frame(Patient, icd)
  nobs <- nrow(d)
  d$code <- rep(NA, nobs)
  
  ## Charlson groups
  
  icd10 <- list()
  
  ################
  ## ICD9 codes ##
  ################
#   ## Myocardial Infarct
#   icd10$codes$mi <- c(410, 412)
#   ## Congestive Heart Failure
#   icd10$codes$chf <- c(428)
#   ## Peripheral Vascular Disease
#   icd10$codes$pvd <- c(4439, 7854)
#   ## Cerebrovascular Disease
#   icd10$codes$cvd <- c(430:438)
#   ## Dementia
#   icd10$codes$dem <- c(290)
#   ## Chronic Pulmonary Disease
#   icd10$codes$cpd <- c(490:496, 500, 501)
#   ## Rheumatic disease
#   icd10$codes$rhm <- c(7100, 7101, 7104, 7140:7142)
#   ## Peptic ulcer disease
#   icd10$codes$pep <- c(531:534)
#   ## Mild liver disease 
#   icd10$codes$mld <- c(5712, 5714:5716)
#   ## Diabetes without chronic complication
#   icd10$codes$dnc <- c(2500:2503, 2507)
#   ## Diabetes with chronic complication
#   icd10$codes$dwc <- c(2504:2506)
#   ## Hemiplegia or paraplegia
#   icd10$codes$ple <- c(342)
#   ## Renal disease 
#   icd10$codes$ren <- c(5830:5832, 5834, 5836, 5837)
#   ## cancer?
#   icd10$codes$can <- c(172)
#   ## Moderate to Severe Liver Disease
#   icd10$codes$liv <- c(5722:5724, 5728)
#   ## Metastatic Solid Tumor
#   icd10$codes$met <- c(196:199)
#   ## AIDS
#   icd10$codes$hiv <- c(42:44)
#   
#   ## ICD9 weights _old_
#   ## from original charlson function
# 
#   icd10$weights <- split(c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,6,6), names(icd10$codes))
  
  ## expand 3 digit ICD10 codes to 4 digits
  ICD3to4 <- function(x) union(unlist(sapply(x, function(y) if(nchar(y)==3){paste(y,0:9, sep="")})), x)
  
  #################
  ## ICD10 codes ##
  #################
  ## Myocardial Infarct
  icd10$codes$mi <- c("I21", "I22", "I252")
  ## Congestive Heart Failure
  icd10$codes$chf <- c("I50")
  ## Peripheral Vascular Disease
  icd10$codes$pvd <- c("I71", "I790", "I739", "R02", "Z958", "Z959")
  ## Cerebrovascular Disease
  icd10$codes$cvd <- c("I60", "I61", "I62", "I63", "I65", "I66", "G450", "G451", "G452", "G458", "G459", "G46",
                 "I64", "G454", "I670", "I671", "I672", "I674", "I675", "I676", "I677", "I678", "I679", "I681", "I682", "I688", "I69")
  ## Dementia
  icd10$codes$dem <- c("F00", "F01", "F02", "F051")
  ## Pulmonary Disease
  icd10$codes$cpd <- c("J40", "J41", "J42", "J44", "J43", "J45", "J46", "J47", "J67", "J44", "J60", "J61", "J62", "J63", "J66", "J64", "J65")
  ## Rheumatic disease/Connective tissue disorder
  icd10$codes$rhm <- c("M32", "M34", "M332", "M053", "M058", "M059", "M060", "M063", "M069", "M050", "M052", "M051", "M353")
  ## Peptic ulcer disease
  icd10$codes$pep <- c("K25", "K26", "K27", "K28")
  ## liver disease 
  icd10$codes$mld <- c("K702", "K703", "K73", "K717", "K740", "K742", "K746", "K743", "K744", "K745")
  ## Diabetes
  icd10$codes$dnc <- c("E109", "E119", "E139", "E149", "E101", "E111", "E131", "E141", "E105", "E115", "E135", "E145")
  ## Diabetes with complications
  icd10$codes$dwc <- c("E102", "E112", "E132", "E142", "E103", "E113", "E133", "E143", "E104", "E114", "E134", "E144")
  ## Paraplegia
  icd10$codes$ple <- c("G81", "G041", "G820", "G821", "G822")
  ## Renal disease 
  icd10$codes$ren <- c("N03", "N052", "N053", "N054", "N055", "N056", "N072", "N073", "N074", "N01", "N18", "N19", "N25")
  ## Cancer
  icd10$codes$can <- c("C0", "C1", "C2", "C3", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C5", "C6", "C70", "C71", "C72", "C73", "C74", "C75", "C76", "C80", "C81", "C82",
                 "C83", "C84", "C85", "C883", "C887", "C889", "C900", "C901", "C91", "C92", "C93", "C940", "C941", "C942", "C943", "C9451", "C947", "C95", "C96")
  ## Metastatic cancer
  icd10$codes$met <- c("C77", "C78", "C79", "C80")
  ## Severe liver disease
  icd10$codes$liv <- c("K729", "K766", "K767", "K721")
  ## AIDS
  icd10$codes$hiv <- c("B20", "B21", "B22", "B23", "B24")
  
  ## ICD10 weights (from V. Sundararajan et al. / Journal of Clinical Epidemiology 57 (2004) 1288–1294)
  
  icd10$weights <- split(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,6), names(icd10$codes))
  
#   ##########################  
#   ## ICD10 codes by score ##
#   ##########################  
#   ## Charlson index (table 1 - Sundararajan et.al. J of Clin. Epi 2004)
#   icd10$codes$Charlson_1point <- c("I21","I22","I252","I50","I71","I790","I739","R02","Z958","Z959","I60","I61","I62","I63","I65","I66","G450","G451","G452","G458","G459","G46","I64","G454","I670","I671","I672","I674","I675","I676","I677","I678","I679","I681","I682","I688","I69","F00","F01","F02","F051","J40","J41","J42","J44","J43","J45","J46","J47","J67","J60","J61","J62","J63","J66","J64","J65","M32","M34","M332","M053","M058","M059","M060","M063","M069","M050","M052","M051","M353","K25","K26","K27","K28","K702","K703","K73","K717","K740","K742","K746","K743","K744","K745","E109","E119","E139","E149","E101","E111","E131","E141","E105","E115","E135","E145")
#   icd10$codes$Charlson_2point <- c("E102","E112","E132","E142","E103","E113","E133","E143","E104","E114","E134","E144","G81","G041","G820","G821","G822","N03","N052","N053","N054","N055","N056","N072","N073","N074","N01","N18","N19","N25","C0","C1","C2","C3","C40","C41","C43","C45","C46","C47","C48","C49","C5","C6","C70","C71","C72","C73","C74","C75","C76","C80","C81","C82","C83","C84","C85","C883","C887","C889","C900","C901","C91","C92","C93","C940","C941","C942","C943","C9451","C947","C95","C96")
#   icd10$codes$Charlson_3point <- c("C77","C78","C79","C80","K729","K766","K767","K721")
#   icd10$codes$Charlson_6point <- c("B20","B21","B22","B23","B24")
#   
#   icd10$codes$Charlson_alcohol <- c("F10","G132","G621","G721","I426","K292","K860","K70","R780","T51","Z721")
#   icd10$codes$HIV_codes <- c("B20","B21","B22","B23","B24")
#   icd10$codes$Pneumococcal_codes <- c("A403","B953", "G001","J13","M001")
  
#   ## ICD10 by score weights
#   
#   icd10$weights <- split(c(1,2,3,6,1,1,1), names(icd10$codes))
  
  
  icd10$codes <- lapply(icd10$codes, ICD3to4)
  
  for (i in 1:nobs) {
    
    patientGroup <- unlist(lapply(icd10$codes, function(x) d[i, "icd"] %in% x))
    d$code[i] <- ifelse(sum(patientGroup)==0, NA, names(icd10$codes)[patientGroup])
  }
  
  mydata <- cast(d, value="icd", Patient ~ code, fun.aggregate=length)
  
  ## omit Patient and NA columns
  chCols <- mydata[, -c(1,ncol(mydata))]>0
  
  comorbidity.n <- rowSums(chCols)
  
  mydata$charlson <- apply((1*chCols), 1, function(x) x%*%unlist(icd10$weights)[colnames(chCols)])
  
  mydata$comorbidity.n <- comorbidity.n
  
  return(mydata)
}

