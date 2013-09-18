#########################################################################################
#
# Functions to extract summary and statistics of interest from the Cox regression output
# then format them into tables for exporting to LaTex and Excel
# and use in reports
#
# Nathan Green
# 11-2012
#
#########################################################################################

extractCox <- function(cox){
  ##
  ## extract a subset of the Cox PH output values
  ## infection status variable assumed to be the last covariate in list
  ## cox: summary(coxph(Surv(start, stop, status)~age+inf, data))
  
  dpl <- 3
  lastRow <- nrow(coef(cox))
  beta <- round(coef(cox)[lastRow,"coef"], dpl)       # exponent of hazard
  se   <- round(coef(cox)[lastRow,"se(coef)"], dpl)       # standard error of beta
  p    <- round(coef(cox)[lastRow,"Pr(>|z|)"], dpl)       # p-value
  CI   <- round(cox$conf.int[lastRow,c("lower .95","upper .95")], dpl)  # lower and upper 95% confidence interval
  
  res <- cbind(beta, "exp(beta)"=round(exp(beta),3), CI[1], CI[2], p)
  
res
}



table.HR <- function(output.HR){
  #
  # Each organism group & Cox PH method alternative format of output data 
  # call: res <- table.HR(output.HR)
  #
  
  namesOrganisms <- names(output.HR)
  namesMethods <- names(output.HR[[1]])
  
  
  colNames <- c("organism", "method", "type", "beta", "exp(beta)", "Lower CI", "Upper CI", "p")
  table.HR <- data.frame(matrix(ncol = length(colNames)))
  
  for (org in namesOrganisms){
    for (meth in namesMethods){
      namesEvent <- names(output.HR[[org]][[meth]]) # different length for different methods
      for (event in namesEvent){
        table.HR <- rbind(table.HR,
                          c(org, meth, event,
                            extractCox(output.HR[[org]][[meth]][[event]])))
      }
    }
  }
  
  colnames(table.HR) <- colNames
  table.HR <- table.HR[!is.na(table.HR[,1]),]		# remove empty rows
  
table.HR
}

#write.table(res, "HCAItable_output.txt", sep="\t")

## Print results in a LaTeX-ready form
#xtable(res)


table2.HR <- function(res, model){
  ##
  ## rearrange table.HR in a report style
  ## used in boc plotting HRboxplot.batch()
  ## model: subdistn, cause-specific
  ##
  ## hr=exp(beta) & upper CI & lower CI
  ## disch time only, disch full, death time only, death full
  
  namesGroup <- unique(res$organism)
  numGroup <- length(namesGroup)
  
  res.sub <- res[res$method==model, c("organism","exp(beta)","Lower CI","Upper CI")]
  
  subHeads <- c("HR","LCI","UCI")
  colHeads <- c("organism", paste("atime",subHeads), paste("afull",subHeads), paste("dtime",subHeads), paste("dfull",subHeads))
  
  res.new <- data.frame(matrix(ncol = length(colHeads), nrow = numGroup), check.rows=FALSE)
  names(res.new) <- colHeads
  
  for (j in 1:numGroup){
    
    res.temp <- NA
    firstrow <- min(which(res.sub$organism==namesGroup[j]))
    
    for (i in 1:4){
      res.temp <- cbind(res.temp, res.sub[firstrow+i-1,-1])
    }
    
    res.new[j,] <- res.temp
  }
  res.new[,1] <- namesGroup
  
res.new
}
## FUNCTION END ##


table3.HR <- function(res, hrtype){
  ## table format used in the JPIDS paper
  ##
  ## organism: is the grouping by organism type or something else
  ## call: table3.HR(res, hrtype="naive")
  ## table3.HR(res, hrtype="timedeptsubdistn")
  ## table3.HR(res, hrtype="timedependentcausespec")
  
  if(hrtype=="timedeptsubdistn"){
    colNames <- c("Group","Disch Time-adjusted","Disch Fully adjusted","Death Time-adjusted","Death Fully adjusted")
  }else if(hrtype=="naive"){
    colNames <- c("Group","Disch", "Death", "Both")
  }else if(hrtype=="timedeptcausespec"){
    colNames <- c("Group","Disch Time-adjusted","Disch Fully adjusted","Death Time-adjusted","Death Fully adjusted","Both Time-adjusted","Both Fully adjusted")
  }
  else {stop("Model type unidentified")}
  
  res.new <- data.frame(matrix(ncol=length(colNames)))
  colnames(res.new) <- colNames
  groupnames <- unique(res$group)
  
  for (name in groupnames){
    
    ## find rows for given organism type and HR method
    whichrows <- which(res$group==name & res$method==hrtype) # & !is.na(res[,"exp(beta)"]))
    rowTotal <- NULL
    for (j in 1:length(whichrows)){
      temp <- paste(res[whichrows[j],"exp(beta)"], " (", res[whichrows[j],"Lower CI"], ", ", res[whichrows[j],"Upper CI"], ")", sep="")
      rowTotal <- c(rowTotal, temp)
    }
    res.new <- rbind(res.new, c(name,rowTotal))
  }
  
  res.new <- res.new[!is.na(res.new[,1]),-1]  	# remove empty rows
  
  if (organism==TRUE){
    ## when group by organism
    ## reformat names and reorder
    rownames(res.new) <- c("All",
                           "Gram-positive",
                           "Gram-negative",
                           "CoNS",
                           "Enterococcus spp.",   
                           "S. aureus",
                           "Other (Gram-positive)",
                           "Other (Gram-negative)",
                           "E. Coli",
                           "Non-p. Streptococci",
                           "Klebsiella spp.",
                           "Enterobacter spp.")
    res.new <- res.new[c(1,2,3,4,9,12,5,11,10,6,7,8),]
  }else{rownames(res.new) <- groupnames}
  
  res.new
}
## END FUNCTION ##


###########################
# Data set summary tables #
###########################

## get interquartile range
iqr <- function(x){paste("[",round(summary(x)[["1st Qu."]],2),",",round(summary(x)[["3rd Qu."]],2),"]", sep="")}

## summary statistics excluding NAs
na.sd <- function(x){sd(x, na.rm=TRUE)}
na.mean <- function(x){mean(x, na.rm=TRUE)}
na.median <- function(x){median(x, na.rm=TRUE)}


summaryTableAll <- function(survData){
  #
  # summary table of dataset descriptive stats by inf/non-inf
  # mean/median (sd or IQR)
  # out <- summaryTableAll(survDataByGroup[["all"]])
  
  survData.mix <- survData[survData$infstatus==0,]  # non-infected patients (controls) only
  survData.inf <- survData[survData$infstatus==1,]  # infected patients (cases) only
  
  dp <- 2
  
  out <- rbind(
    
    ## sample sizes
    c(format(nrow(survData.inf),nsmall=1), round(nrow(survData.inf)/nrow(survData),dp),
      round(nrow(survData.mix)), round(nrow(survData.mix)/nrow(survData),dp),
      round(nrow(survData)), 1),
    
    ## ages
    #c(paste(round(mean(survData.inf$age)),"/",median(survData.inf$age)),round(sd(survData.inf$age)),
    c(paste(round(mean(survData.inf$age)),"/",median(survData.inf$age)),iqr(survData.inf$age),
      #paste(round(mean(survData.mix$age)),"/",median(survData.mix$age)),round(sd(survData.mix$age)),
      paste(round(mean(survData.mix$age)),"/",median(survData.mix$age)),iqr(survData.mix$age),
      #paste(round(mean(survData$age)),"/",median(survData$age)),round(sd(survData$age))),
      paste(round(mean(survData$age)),"/",median(survData$age)), iqr(survData$age) ),
    
    ## length of stays
    #c(paste(round(mean(survData.inf$time)),"/",median(survData.inf$time)),round(sd(survData.inf$time)),
    #  paste(round(mean(survData.mix$time, na.rm=T)),"/",median(survData.mix$time, na.rm=T)),round(sd(survData.mix$time, na.rm=T)),
    #  paste(round(mean(survData$time, na.rm=T)),"/",median(survData$time, na.rm=T)),round(sd(survData$time, na.rm=T))),
    c(paste(round(mean(survData.inf$time),dp),"/",median(survData.inf$time)),iqr(survData.inf$time),
      paste(round(mean(survData.mix$time, na.rm=T),dp),"/",median(survData.mix$time, na.rm=T)),iqr(survData.mix$time),
      paste(round(mean(survData$time, na.rm=T),dp),"/",median(survData$time, na.rm=T)),iqr(survData$time) ),
    
    ## infection times
    #c(paste(round(mean(survData.inf$spectime)),"/",median(survData.inf$spectime)),round(sd(survData.inf$spectime)),
    # paste(round(mean(survData.mix$spectime)),"/",median(survData.mix$spectime)),round(sd(survData.mix$spectime)),
    #paste(round(mean(survData$spectime)),"/",median(survData$spectime)),round(sd(survData$spectime))),  
    c(paste(round(mean(survData.inf$spectime),dp),"/",median(survData.inf$spectime)),iqr(survData.inf$spectime),
      paste(round(mean(survData.mix$spectime),dp),"/",median(survData.mix$spectime)),iqr(survData.mix$spectime),
      paste(round(mean(survData$spectime),dp),"/",median(survData$spectime)),iqr(survData$spectime) ), 
    
    ## in-hospital deaths
    c(round(table(survData.inf$event)[2]),round(table(survData.inf$event)[2]/(table(survData.inf$event)[2]+table(survData.inf$event)[1]),dp),
      round(table(survData.mix$event)[2]),round(table(survData.mix$event)[2]/(table(survData.mix$event)[2]+table(survData.mix$event)[1]),dp),
      round(table(survData$event)[2]),round(table(survData$event)[2]/(table(survData$event)[2]+table(survData$event)[1]),dp)),
    
    ## sex (female)
    c(round(table(survData.inf$gender)[1]),round(table(survData.inf$gender)[1]/(table(survData.inf$gender)[1]+table(survData.inf$gender)[2]),dp),
      round(table(survData.mix$gender)[1]),round(table(survData.mix$gender)[1]/(table(survData.mix$gender)[1]+table(survData.mix$gender)[2]),dp),
      round(table(survData$gender)[1]),round(table(survData$gender)[1]/(table(survData$gender)[1]+table(survData$gender)[2]),dp)))
  
  rownames(out) <- c("Patient sample size","Age (years)","LoS (days)","Time from admission to infection (days)","Deaths (frequency)","Sex (F) (frequency)")
  colnames(out) <- c("HA-BSI", "Prop", "Non-HA-BSI", "Prop", "All", "Prop")
  
  ## rearrange columns
  out <- out[,c(3,4,1,2,5,6)]
  
  # write.table(out, ".\\output\\summaryTableAll.txt")
  
  return(pandoc.table(out, caption="Table: Dataset summary statistics including risk factors, comorbidites and patient movements summary statistics.
                      For count dat the subset size if given; for continuous values mean/median is given.
                      Note that a patient can be in more than one risk factor group.", style = "grid", split.tables=Inf, justify="left"))
}
## END FUNCTION ##



summaryTableGroup <- function(survDataByGroup){
  ##
  ## summary table of descriptive statistics by (organism) group
  ## sd or IQR
  ## call: out <- summaryTableOrg(survDataByGroup)
  
  out <- NA
  dp <- 2
  
  for (group in names(survDataByGroup)){
    out <- rbind(out,
                 c(group,
                   ##sample size
                   nrow(survDataByGroup[[group]][survDataByGroup[[group]]$infstatus==1,]),
                   
                   ## age
                   paste(round(na.mean(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus==1]),dp),"/",
                         na.median(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus==1]),
                         #   " (",round(na.sd(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus==1])),")",sep=""),
                         iqr(survDataByGroup[[group]]$age[survDataByGroup[[group]]$infstatus==1]), sep=""),
                   
                   ## gender
                   round(table(survDataByGroup[[group]]$gender[survDataByGroup[[group]]$infstatus==1])[1],dp),
                   
                   ## LoS
                   paste(round(na.mean(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus==1]),dp),"/",
                         na.median(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus==1]),
                         # " (",round(na.sd(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus==1])),")", sep=""),
                         iqr(survDataByGroup[[group]]$time[survDataByGroup[[group]]$infstatus==1]), sep=""),
                   
                   ## infection time
                   paste(round(na.mean(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus==1]),dp),"/",
                         na.median(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus==1]),
                         # " (",round(na.sd(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus==1])),")",sep=""),
                         iqr(survDataByGroup[[group]]$spectime[survDataByGroup[[group]]$infstatus==1]), sep=""),
                   
                   ## deaths
                   round(table(survDataByGroup[[group]]$event[survDataByGroup[[group]]$infstatus==1])[2]),
                   
                   sum(survDataByGroup[[group]]$cancer[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$prem[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$cong[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$surgical[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$cath[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$Tai[survDataByGroup[[group]]$infstatus==1]),
                   sum(survDataByGroup[[group]]$highRisk[survDataByGroup[[group]]$infstatus==1])
                 ))
    
  }
  colnames(out) <- c("Organism","Sample size","Age (years)","Sex (F)","LoS (days)","Time from admission to infection (days)","Deaths",
                     "Cancer","Premature birth","Congenital disease","Surgical","In-dwelling catheter","Tai","At least one risk factor")
  rownames(out) <- out[,"Organism"]
  out <- out[!is.na(out[,1]),-1]    # remove empty rows
  
  ## rearrange rows
  #out <- out[c("all", "1", "-1", "COAGULASE NEGATIVE STAPHYLOCOCCUS", "E. COLI", "ENTEROBACTER", "ENTEROCOCCUS", "KLEBSIELLA", "NON-PYOGENIC STREPTOCOCCUS",
  #             "STAPHYLOCOCCUS AUREUS", "other",  "P. AERUGINOSA", "MICROCOCCUS", "STREP B", "SALMONELLA", "STREPTOCOCCUS PNEUMONIAE", "N. MENINGITIDIS", "STREP A", "ACINETOBACTER"),]                                                 
  out <- out[c("all", "1", "-1", "COAGULASE NEGATIVE STAPHYLOCOCCUS", "E. COLI", "ENTEROBACTER", "ENTEROCOCCUS", "KLEBSIELLA", "NON-PYOGENIC STREPTOCOCCUS",
               "STAPHYLOCOCCUS AUREUS", "other (Gram-positive)", "other (Gram-negative)",  "P. AERUGINOSA", "MICROCOCCUS", "STREP B", "SALMONELLA", "STREPTOCOCCUS PNEUMONIAE", "N. MENINGITIDIS", "STREP A", "ACINETOBACTER"),]                                                 
  
  #write.table(out, ".\\output\\summaryTableOrg.txt")
  
  return(out)
  #  pandoc.table(out, caption="caption:...", style = "grid")
}
## END FUNCTION ##



summaryTableRF <- function(survData){
  ##
  ## output summary table of the patients
  ## comorbidities and risk factors
  ## split by infected and non-infected cases
  ##
  ## call: out <- summaryTableRF(survDataByGroup$all)
  
  out <- NULL
  emptyRow <- c(NA, NA, NA, NA, NA, NA)
  dp <- 2
  
  ncase <- nrow(survData)
  
  ## risk factor and comorbidities
  ## record true (present) cases only
  
  ## empty rows for extra row labels
  out <- rbind(out, emptyRow)
  
  freqRow <- function(survData, out, rf){
    x <- as.data.frame(table(survData$infstatus, survData[,rf]))
    y <- as.data.frame(prop.table(table(survData$infstatus, survData[,rf]),1))
    all <- c(sum(survData[,rf]), round(sum(survData[,rf])/length(survData[,rf]), dp)) 
    out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
    out
  }
  
  out <- freqRow(survData, out, "cancer")
  out <- freqRow(survData, out, "prem")
  out <- freqRow(survData, out, "cong")
  out <- freqRow(survData, out, "surgical")
  out <- freqRow(survData, out, "cath")
  out <- freqRow(survData, out, "Tai")
  out <- freqRow(survData, out, "highRisk")
  
  
  ## type of admission (admission method)
  out <- rbind(out, emptyRow)
  
  ### Elective
  elective <- c("Elective - booked", "Elective - planned", "Elective - from waiting list")
  x <- as.data.frame(with(survData,
                          table(infstatus, hes_admimethdescription%in%elective)))
  y <- as.data.frame(prop.table(with(survData,
                                     table(infstatus, hes_admimethdescription%in%elective)),1))
  all <- c(sum(survData$hes_admimethdescription%in%elective), round(sum(survData$hes_admimethdescription%in%elective)/ncase, dp))
  out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
  
  
  ### Emergency
  emergency <- c("Emergency - other means, including patients who arrive via A&E department of another HC provider",
                 "Emergency - via A&E services, including casualty department of provider",
                 "Emergency - via General Practitioner (GP)",
                 "Emergency - via Bed Bureau, including Central Bureau",
                 "Emergency - via consultant out-patient clinic")
  
  x <- as.data.frame(with(survData,
                          table(infstatus, hes_admimethdescription%in%emergency)))
  y <- as.data.frame(prop.table(with(survData,
                                     table(infstatus, hes_admimethdescription%in%emergency)),1))
  all <- c(sum(survData$hes_admimethdescription%in%emergency), round(sum(survData$hes_admimethdescription%in%emergency)/ncase, dp))
  out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
  
  ## total
  out <- rbind(out, aggregate(as.numeric(tail(out,2)), by=list(c(1,1,2,2,3,3,4,4,5,5,6,6)), sum)$x)
  
  ## intensive neocare
  neocare <- c("Level 1 intensive care", "Level 2 intensive care")
  x <- as.data.frame(with(survData,
                          table(infstatus, hes_neocaredescription%in%neocare)))
  y <- as.data.frame(prop.table(with(survData,
                                     table(infstatus, hes_neocaredescription%in%neocare)),1))
  all <- c(sum(survData$hes_neocaredescription%in%neocare), round(sum(survData$hes_neocaredescription%in%neocare)/ncase, dp))
  out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
  
  
  
  ## origin of patient (admission description)
  out <- rbind(out, emptyRow)
  
  ### another hospital
  transfer.txt <- "Transfer of any admitted patient from another hospital provider"
  x <- as.data.frame(with(survData,
                          table(infstatus, hes_admimethdescription==transfer.txt)))
  y <- as.data.frame(prop.table(with(survData,
                                     table(infstatus, hes_admimethdescription==transfer.txt)),1))
  all <- c(sum(survData$hes_admimethdescription==transfer.txt), round(sum(survData$hes_admimethdescription==transfer.txt)/ncase, dp))
  out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
  
  ### residence
  residence <- c("The usual place of residence, including no fixed abode",
                 "Temporary place of residence when usually resident elsewhere")
  x <- as.data.frame(with(survData,
                          table(infstatus, hes_admisorcdescription%in%residence)))
  y <- as.data.frame(prop.table(with(survData,
                                     table(infstatus, hes_admisorcdescription%in%residence)),1))
  all <- c(sum(survData$hes_admisorcdescription%in%residence), round(sum(survData$hes_admisorcdescription%in%residence)/ncase, dp))
  out <- rbind(out, c(x[x$Var2==T,"Freq"], round(y[y$Var2==T,"Freq"],dp), all))
  
  ## total
  out <- rbind(out, aggregate(as.numeric(tail(out,2)), by=list(c(1,1,2,2,3,3,4,4,5,5,6,6)), sum)$x)
  
  ## rearange columns
  if (ncol(out)==6){
    colnames(out) <- c("Non-HA-BSI", "HA-BSI", "Prop", "Prop", "All", "Prop")
    out <- out[,c(1,3,2,4,5,6)]} 
  else {colnames(out) <- c("Non-HA-BSI", "Prop")}
  
  rownames(out) <- c("Risk Factors", "Cancer", "Premature birth", "Congenital disorder", "Surgical", "In-dwelling catheter", "Tai", "At least one",
                     "Type of Admission", "Elective", "Emergency", "TOTAL",
                     "Intensive neonatal care",
                     "Origin of patient", "Another hospital", "Residence", "TOTAL")
  
  
  #write.table(out, ".\\output\\summaryTableRF.txt")
  
return(out)
#return(pandoc.table(out, caption="caption:...", style = "grid"))
}
## END FUNCTION ##


LOStable <- function(output.LOS, se=TRUE, orgs=TRUE){
  ##
  ## excess length of stay table by group
  ## se: standard error or confidence interval
  ## orgs: is grouping by organsim type or not
  ## call: LOS <- LOStable(output.LOS)
  
  LOS <- NULL
  dp <- 2
  
  for (i in names(output.LOS)){
    if(is.na(output.LOS[[i]][[1]][1])){ LOS <- c(LOS, NA)
    } else {
      if(se){
        LOS <- c(LOS, paste(round(output.LOS[[i]]$clos.data$e.phi,dp), " (", round(output.LOS[[i]]$se,dp), ")", sep="")  )
      }else{
        LOS <- c(LOS, paste(round(output.LOS[[i]]$clos.data$e.phi,dp), " (",
                            round(output.LOS[[i]]$clos.data$e.phi-1.96*output.LOS[[i]]$se,dp), ", ",
                            round(output.LOS[[i]]$clos.data$e.phi+1.96*output.LOS[[i]]$se,dp), ")", sep="")  )
      }
    }
  }
  
  #res <- data.frame(Group=capwords(names(output.LOS), strict = TRUE), LOS)
  res <- data.frame(Group=names(output.LOS), "Excess LOS"=LOS)
  
  if(orgs==TRUE){
    ## for specific organism grouping
    ## format names and reorder
    rownames(res) <- c("All",
                       "Gram-positive",
                       "Gram-negative",
                       "CoNS",
                       "Enterococcus spp.",   
                       "S. aureus",
                       "Other (Gram-positive)",
                       "Other (Gram-negative)",
                       "E. Coli",
                       "Non-p. Streptococci",
                       "Klebsiella spp.",
                       "Enterobacter spp.")
    neworder <- c(1,2,3,4,9,12,5,11,10,6,7,8)
  }else{
    ## general, generic default grouping
    rownames(res) <- names(output.LOS)
    neworder <- 1:NROW(res)
  }    
  
  res <- res[neworder, -1, drop=FALSE]
  
  #write.table(res, ".\\output\\tables\\LOStable.txt", sep="\t")
res
}
## END FUNCTION ##


commonlist <- function(group, cut=20, plot=FALSE){
  #
  # Function that produces an ordered list by frequency.
  # cut: frequency cut-off value
  #
  
  x <- as.data.frame(table(group))  		# table of frequencies
  x <- data.frame(x, prop=x$Freq/sum(x$Freq))	# percentage of total	
  
  if (plot){
    plot(x$Var1[x$Freq>cut], x$Freq[x$Freq>cut], las=3)}
  
  x[x$Freq>cut,] [order(x[x$Freq>cut,"Freq"], decreasing=TRUE),]	# cutoff and reorder
}
## END FUNCTION ##


HCAIsummary <- function(survData){
  ##
  ## Function to aggregate the operations for summary statistics for age, length of stay & BSI time
  ## output to console
  
  trimPerc <- 0.0        # trimmed/Winsorising mean cut-off percentage
  
  print(c("age", summary(survData$age)))					# quartiles
  print(sd(survData$age, na.rm=TRUE)) 
  print(table(survData$age))           # frequencies  
  
  print(c("LOS", summary(survData$time)))        # hospital length of stay
  print(c("sd", sd(survData$time, na.rm=TRUE)))
  print(c("trimmed mean", mean(survData$time, trim=trimPerc)))			# trimmmed mean
  #winsor.means(survData$time, trim = trimPerc)  # winsorised mean comparison
  
  print(c("spectime", summary(survData$spectime)))
  print(c("sd", sd(survData$spectime, na.rm=TRUE)))
  
  print(c("deaths", table(survData$event)))   # proportion of death-discharge
  print(c("gender", table(survData$gender)))	# proportion of male-female
  
return()
}
## END FUNCTION ##



catcodes <- function(codes){
  #
  # rearrange the reference look-up table
  # so that codes in a single column
  #
  # codes <- read.csv(".\\Reference_tables\\explicitSMEcodes.csv")
  # codes <- read.csv(".\\Reference_tables\\SSI_OPCS_codes.csv")
  # codes <- read.table(".\\Reference_tables\\ICDGroup_match.txt")

  codes$Code <- clean.ref(codes$Code) 

  x <- aggregate(codes, by=list(codes$Group), paste, collapse=",")[,1:2]
  write.table(x, file=".\\Reference_tables\\matchcodesGrouped.txt", sep="\t")
}
## END FUNCTION ##