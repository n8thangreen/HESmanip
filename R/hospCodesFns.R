#----------------------------------
#
#	ICD10 and OPCS functions
# cleaning, searching and grouping
#
# Nathan Green
# 12-2012
#
#-----------------------------------



clean.data <- function(data, cols){
#
# in ICD10 codes remove unwanted trailing characters
#
	for (i in cols){
		data[,i] <- gsub("-$","0",data[,i])
		data[,i] <- gsub("X$","0",data[,i])
		data[,i] <- gsub("X.$","0",data[,i])
		data[,i] <- substring(data[,i],1,4)
    
    #data[,i] <- gsub("X.$","0",data[,i])
		#data[,i] <- gsub("O03.","O03",data[,i])
		#data[,i] <- gsub("N04.","N04",data[,i])
		#data[,i] <- gsub("N05.","N05",data[,i])
	}
data
}
## END FUNCTION ##


clean.ref <- function(refcodes){
#
# clean reference table codes
#
	refcodes <- as.character(refcodes)
	refcodes <- gsub("'","",refcodes)
	refcodes <- gsub(" ","",refcodes)
	refcodes <- gsub(" ","",refcodes)
	refcodes <- gsub(" ","",refcodes)
	refcodes <- gsub("+","",refcodes, fixed=TRUE)
	refcodes <- gsub("*","",refcodes, fixed=TRUE)
	refcodes <- gsub(".","",refcodes, fixed=TRUE)
	refcodes <- gsub("-","",refcodes, fixed=TRUE)

refcodes
}
## END FUNCTION ##


ICDwordmatch <- function(ICD10){
  ##
  ## create a list of ICD10 codes associated with key words
  ## ICD10: lookup table
  ## RFcodes: lists of codes for each key word
  
  RFcodes <- list()
  #cat <- c("cancer", "cathet", "dialys", "low birth", "preterm", "prematurity", "^intravenous")
  cat <- c("cancer","cathet","dialys","low birth","preterm","prematurity","^intravenous",
           " sepsis"," septic"," infection")   #"bacteremia"  ## with additional infection key words (include leading space to avoid eg "aseptic")
  
  for (i in seq_along(cat)){
    ## logical vector of all patients if match in either "Group" or "Description" field
    
    catmod <- paste("(?i)",cat[i],sep="")   # 'i' ignores lower or upper case (insensitive)
    RFrows <- unique(
      c(grep(catmod, ICD10$Group, fixed=FALSE),
        grep(catmod, ICD10$Description, fixed=FALSE)))  # not exact string matching permitted
    
    RFcodes[[i]] <- ICD10$Code[RFrows]  # vector of codes for each category in the list
  }
  
  names(RFcodes) <- cat     # list category names
  
RFcodes
}
## END FUNCTION ##



whichPatientDrug <- function(inf.data, drugGuideline.data, patientdrugs){
  ## Determine which patients are susceptible or resistant to the guideline drugs,
  ## accounting for age category
  ##
  ## patientdrugs: list of drugs per patient
  ## returns: list of susc/res status per patient a) per drug b) over all drugs
  ##          TRUE entry if guideline antibiotic and resistant/susceptible
  ##
  ## call: whichPatientDrug(patientdrugs.res)
  
  npatient <- length(patientdrugs)            # number of patients
  ndrugs <- length(drug.names)                # number of antibiotics (names from LabBase2 column headings)
  nguidelines <- nrow(drugGuideline.data)     # number of guidelines
  
  drug.names.orig <- drug.names
  drug.names <- toupper(drug.names)
  
  FollowGuideline.tot   <- rep(FALSE, npatient)                   # all drugs combined for single result per patient
  FollowGuideline.indiv <- matrix(NA, nrow=npatient, ncol=ndrugs) # individual drug results per patient
  
  colnames(FollowGuideline.indiv) <- drug.names
  
  ## place patients in corresponding age category to look-up table
  patient.age_cat <- AgeCatModify(inf.data)
  
  
  for (i in 1:nguidelines){
    
    ## check if guideline drug field is not empty and relevant to dataset
    ## check against individual and generic (grouped) codes
    if (drugGuideline.data$Antimicrobial.Code[i]!="" & 
          (toupper(drugGuideline.data$Antimicrobial.Code[i])%in%drug.names |
             toupper(drugGuideline.data$Generic.Antimicrobial.Code)[i]%in%drug.names)){
      
      ## NB duplicate entries will appear for same antibiotics for different ages
      
      ## whether ageCat is included in look-up table
      ageMatch <- ifelse("ageCat"%in%names(drugGuideline.data), patient.age_cat%in%c(drugGuideline.data$ageCat[i],""), TRUE)
      
      ## for each patient TRUE/FALSE if match in guidelines fields
      # 1 x npatient vector
      drugMatch <- sapply(patientdrugs, function(drug.list)
        ((toupper(drugGuideline.data$Antimicrobial.Code[i])%in%toupper(drug.list)) |
           (toupper(drugGuideline.data$Generic.Antimicrobial.Code[i])%in%toupper(drug.list)))) &
        ageMatch &
        ## need to ensure that entries in different tables are consistent in order to match pathogen names
        ## case-insensitive matching
        (toupper(drugGuideline.data$suspected_organism[i]) == toupper(inf.data$lab_organismname) |
           toupper(drugGuideline.data$suspected_organism[i]) == toupper(inf.data$lab_group) |
           drugGuideline.data$suspected_organism[i] == "")
      ## use logical OR because same drug may appear in guideline more than once
      ## flag antibiotic if TRUE
      drugCol <- ifelse(toupper(drugGuideline.data$Antimicrobial.Code[i])%in%drug.names,
                        toupper(as.character(drugGuideline.data$Antimicrobial.Code[i])), toupper(as.character(drugGuideline.data$Generic.Antimicrobial.Code[i])))
      FollowGuideline.indiv[,drugCol] <- FollowGuideline.indiv[,drugCol] | drugMatch 
      FollowGuideline.tot <- FollowGuideline.tot | drugMatch
    }
  }
  
  ## replace fields with original (non-capitalised) names
  drug.names <- drug.names.orig
  colnames(FollowGuideline.indiv) <- drug.names.orig
  
list(indiv=FollowGuideline.indiv, all=FollowGuideline.tot)
}
## END FUNCTION ##


whichPatientDrug.joins <- function(inf.drug.data, drugGuideline.data){
  ## alternative to whichPatientDrug() using (melt &) joins instead of loops
  ## returns different, combined format (list) to whichPatientDrug()
  ## 
  ## inf.drug.data: array of patient by drug names containing 0,1 or NA
  ## drugGuideline.data: lookup table of drug-bug pairs of first-line treatment
  
  ## TODO ##
  ## to finish and test
  
  
  lookup <- drugGuideline.data[,c("suspected_organism","Generic.Antimicrobial.Code")] #"ageCat"
  colnames(lookup) <- c("organism", "drug") #"ageCat"
  lookup <- lookup[!duplicated(lookup),]
  
  ## remove "lab_" prefix
  colnames(inf.drug.data) <- drug.names
  ## patient ids
  inf.drug.data$patient <- rownames(inf.drug.data)
  inf.drug.data <- cbind(inf.drug.data, organism=inf.data$lab_group)
  # convert to long format
  moltendata <- melt(inf.drug.data, id=c("patient", "organism"))
  colnames(moltendata) <- c("patient", "organism", "drug", "value")
  moltendata <- moltendata[order(moltendata$patient),]
  
  ## case-insensitive matching
  lookup <- apply(lookup, 2, toupper); moltendata <- apply(moltendata, 2, toupper)
  
  ## filter by guideline drugs
  ## antibiotics for _named_drugs_
  GLoutput <- merge(moltendata, lookup, by=c("drug", "organism"))
  GLoutput <- subset(GLoutput, select=c("drug","patient","value"))
  
  ## antibiotics for _no_specific_organism_
  GLoutput.generic <- merge(moltendata, lookup[lookup[,"organism"]=="",], by="drug")
  GLoutput.generic <- subset(GLoutput.generic, select=c("drug","patient","value"))
  
  ## combine
  GLoutput <- rbind(GLoutput, GLoutput.generic)
  
  GLoutput$patient <- as.numeric(as.character(GLoutput$patient))
  GLoutput <- GLoutput[order(GLoutput$patient),]
  GLoutput <- GLoutput[!duplicated(GLoutput),]
  
  
  ## TODO ##
  ## join including age category of patients
  
  
  GLstrata <- list()
  GLstrata$res <- GLoutput[GLoutput$value==1 & !is.na(GLoutput$value),]
  GLstrata$sus <- GLoutput[GLoutput$value==0 & !is.na(GLoutput$value),]
  GLstrata$notest <- GLoutput[is.na(GLoutput$value),]
  
return(GLstrata)
}
## END FUNCTION ##


getPatientGroupings <- function(inf.data, mix.data){
  ##
  ## use a combination of fields identified through partial word matching
  ## and fields explicitly identified by subject matter experts and look-up table sources
  ## inf.data: infected patient data
  ## mix.data: non-infected patient data
  
  
  ## codes identified by subject matter experts
  SMEcodes <- read.csv(".\\Reference_tables\\explicitSMEcodes.csv")
  #SMEcodes <- read.csv(".\\Reference_tables\\explicitSMEcodes_temp.csv")
  
  ###########################################################
  # OPCS Classification of Interventions & Procedures codes #
  # group by operation and procedure patient codes OPCS     #
  ###########################################################
  
  codeColnames <- c("Code", "Group", "Description", "Scheme")
  
  ## HPA Surgical Site Infection Surveillance Service Report, used to identify "surgical" patients
  SSI_OPCS_codes <- read.csv(".\\Reference_tables\\SSI_OPCS_codes.csv")
  colnames(SSI_OPCS_codes) <- codeColnames                                    # renames columns to match other tables
  SSI_OPCS_codes$Code <- clean.ref(SSI_OPCS_codes$Code)  
  
  ## OPCS codes columns in main HES/LabBase2 main dataset
  opcode_cols.inf <- grep("hes_op([0-9]+)_(main)?code", names(inf.data))
  opcode_cols.mix <- grep("operativeprocedure", names(mix.data))
  
  inf.data <- clean.data(inf.data, opcode_cols.inf)
  mix.data <- clean.data(mix.data, opcode_cols.mix)
  
  ## surgical
  surgFlags.inf <- facFlags(SSI_OPCS_codes, "OPCS", inf.data[,opcode_cols.inf])
  surgFlags.mix <- facFlags(SSI_OPCS_codes, "OPCS", mix.data[,opcode_cols.mix])
  
  ## bespoke groups
  SMEcodesFlags.inf <- facFlags(SMEcodes, "OPCS", inf.data[,opcode_cols.inf])
  SMEcodesFlags.mix <- facFlags(SMEcodes, "OPCS", mix.data[,opcode_cols.mix])
  
  
  
  #########################
  # ICD10 codes           #
  # group by disease code #
  #########################
  
  diagcode_cols.mix <- seq(13, 51, by=2)  # extract columns of ICD10 code from main HES/LabBase2 database ##TODO## not straightforward because empty column headers as v##
  diagcode_cols.inf <-  grep("hes_diag([0-9]+)_code", names(inf.data))   # match diagnosis text in column names to identify
  
  ## clean ICD10 data: e.g. remove unwanted trailing characters
  inf.data <- clean.data(inf.data, diagcode_cols.inf)
  mix.data <- clean.data(mix.data, diagcode_cols.mix)
  
  ## combine SME output for ICD10 and OPCS
  SMEcodesFlags.inf <- SMEcodesFlags.inf | facFlags(SMEcodes, "ICD10", inf.data[,diagcode_cols.inf])
  SMEcodesFlags.mix <- SMEcodesFlags.mix | facFlags(SMEcodes, "ICD10", mix.data[,diagcode_cols.mix])
  
  
  
  #   
  #   #source(".\\HCAI_R_code\\icd_function_v0.91.R")  # list of ICD10 by subchapter groups, returns: icd_codes_subchap
  #   #source(".\\HCAI_R_code\\icd_codes_chap.R")  	  # list of ICD10 by chapter, returns: icd_codes_chap
  #   ICD10 <- read.csv(".\\Reference_tables\\ccs_icd10.csv", header=TRUE)	# H-CUP US Clinical Classification Software groupings array
  #                                                                         # http://www.hcup-us.ahrq.gov/toolssoftware/icd_10/ccs_icd_10.jsp
  #   colnames(ICD10) <- codeColnames[1:3]
  #   
  #   ## clean ICD10 reference table
  #   ICD10$Code <- clean.ref(ICD10$Code)
  #     
  #   ##### word matching ####
  #   ## select grouping ref 
  #   #ICDGroup <- icd_codes_subchap
  #   #ICDGroup <- icd_codes_chap
  #   #ICDGroup <- ICD10		          	# group the ICD10 codes according to H-CUP Clinical Classifications Software (CCS)
  #   
  #   ## 1) create lookup table from matching
  #   #ICDGroup.match <- ICDwordmatch(ICD10)   # create lists of ICD10 codes for each key word
  #   #ICDGroup.match <- melt.list(ICDGroup.match) # as.array to coerse to a common format
  #   #colnames(ICDGroup.match) <- codeColnames[1:2]
  #   #ICDGroup.match[,2] <- as.factor(ICDGroup.match[,"Group"])
  #   #ICDGroup.match$Scheme <- "ICD10"
  #   
  #   ## 2) or read in preprocessed matching lookup table
  ICDGroup.match <- read.table(".\\Reference_tables\\ICDGroup_match.txt")
  
  ## matrix of indicators for each group and each patient
  ICDgroupFlag.mix <- facFlags(ICDGroup.match, "ICD10", mix.data[,diagcode_cols.mix])
  ICDgroupFlag.inf <- facFlags(ICDGroup.match, "ICD10", inf.data[,diagcode_cols.inf])
  
  
  ## combine OPCS codes with relevant ICD10 fields
  cathCodesgroupFlags.inf <- SMEcodesFlags.inf[,"invasive"] | apply(ICDgroupFlag.inf[,c("cathet","dialys","^intravenous")],1,any)
  cathCodesgroupFlags.mix <- SMEcodesFlags.mix[,"invasive"] | apply(ICDgroupFlag.mix[,c("cathet","dialys","^intravenous")],1,any)
  
  premCodesgroupFlags.inf <- SMEcodesFlags.inf[,"low birth"] | apply(ICDgroupFlag.inf[,c("low birth","preterm","prematurity")],1,any)
  premCodesgroupFlags.mix <- SMEcodesFlags.mix[,"low birth"] | apply(ICDgroupFlag.mix[,c("low birth","preterm","prematurity")],1,any)
  
  infCodesgroupFlags.inf <- SMEcodesFlags.inf[,"infection"] #| apply(ICDgroupFlag.inf[,c(" infection", " sepsis", " septic")],1,any)
  infCodesgroupFlags.mix <- SMEcodesFlags.mix[,"infection"] #| apply(ICDgroupFlag.mix[,c(" infection", " sepsis", " septic")],1,any)
  
  codesIndicators.inf <- data.frame(surg=surgFlags.inf,
                                    cong=SMEcodesFlags.inf[,"congenital"],
                                    head=SMEcodesFlags.inf[,"head"],
                                    infection=infCodesgroupFlags.inf,
                                    cath=cathCodesgroupFlags.inf,
                                    prem=premCodesgroupFlags.inf,
                                    mental=SMEcodesFlags.inf[,"mental"],
                                    other=SMEcodesFlags.inf[,"other"],
                                    Tai=SMEcodesFlags.inf[,"Tai"],
                                    unknown=SMEcodesFlags.inf[,"unknown"],
                                    cancer=ICDgroupFlag.inf[,"cancer"])
  
  codesIndicators.mix <- data.frame(surg=surgFlags.mix,
                                    cong=SMEcodesFlags.mix[,"congenital"],
                                    head=SMEcodesFlags.mix[,"head"],
                                    infection=infCodesgroupFlags.mix,
                                    cath=cathCodesgroupFlags.mix,
                                    prem=premCodesgroupFlags.mix,
                                    mental=SMEcodesFlags.mix[,"mental"],
                                    other=SMEcodesFlags.mix[,"other"],
                                    Tai=SMEcodesFlags.mix[,"Tai"],
                                    unknown=SMEcodesFlags.mix[,"unknown"],
                                    cancer=ICDgroupFlag.mix[,"cancer"])
  
  
  ## combine all covariates into a high risk group
  codesIndicators.inf <- cbind(codesIndicators.inf, highRisk=rowSums(codesIndicators.inf)>0)
  codesIndicators.mix <- cbind(codesIndicators.mix, highRisk=rowSums(codesIndicators.mix)>0)
  
list(inf=codesIndicators.inf, mix=codesIndicators.mix)
}
## FUNCTION END ##


my_recode <- function(codes, lookuplist){
##
## groups against a list
## using change in factor labels

  nfac <- factor(codes)
  lfac <- levels(nfac)
  
  othrlevs <- lfac[ !lfac %in% unlist(lookuplist) ]
  
  levels(nfac)<- c(lookuplist, list(all_others=othrlevs))
  
nfac
}
## END FUNCTION ##


facFlags <- function(lookup, scheme, data){
  ##
  ## group patients using change in factor labels
  ## over all groups
  
  ## convert look-up array to list
  ## need as.list() in case of only one group
  lookup.list <- as.list(unstack(lookup[lookup$Scheme==scheme, c("Code","Group")], Code~Group))
  
  ## recode by group each column of codes
  y <- apply(data, 2, function(x) my_recode(x, lookup.list))
  list.names <- names(lookup.list)
  res <- NULL
  
  ## look across all codes for a given patients (columns) and combine
  for (i in list.names){
    res <- cbind(res, apply(y, 1, function(x) any(x==i)) )
  }
  colnames(res) <- list.names
  
res
}
## END FUNCTION ##


AgeCatModify <- function(inf.data){
  ## find which corresponding age category in the guidelines
  ## are the patients in the LabBase2 database
  ## age categories different in guidelines and databases so
  ## need to regroup & relabel
  
  #############################################
  # 1) <48 hours 7001                         #    
  # 2) 48hrs - 1 month 7002:7003              #
  # 3) 1 month - 18 years c(7004:7007, 1:18)  #
  # 4) over 18 19:200                         #
  #############################################
  
  patient.age_cat <<- cut(inf.data[,"hes_ageatstartofepisode"], breaks=c(0,18,200,7001.5,7003.5,7008))
  patient.age_cat <- factor(patient.age_cat, c("(200,7002]", "(7002,7004]", "(7004,7008]", "(0,18]", "(18,200]")) # reorder
  
  patient.age_cat[patient.age_cat=="(0,18]"] <- "(7004,7008]"   # combine 2 age categories
  patient.age_cat <- droplevels(patient.age_cat)  # or factor(patient.age_cat)
  levels(patient.age_cat) <- 1:4  #levels(drugGuideline.data[,"ageCat"] )  # rename
  
patient.age_cat
}
## END FUNCTION ##



