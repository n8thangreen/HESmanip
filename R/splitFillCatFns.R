splitByOrg <- function(total.data, survData.total){
  ## split data into 2 predefined organism groups for later HR
  ## with particular infection type or not
  ## other infection types removed and do not go in with everything else
  ##
  ## call: survDataByGroup <- splitByOrg(total.data, survData.total)
  ##
  ## DEPRECATED by
  ## c(splitByGroup(total.data, survData.total, groupCol="lab_group"),
  ##   splitByGroup(total.data, survData.total, groupCol="Gram"))
  
  survDataByGroup <- list()
  namesGroup <- unique(total.data$lab_group[!is.na(total.data$lab_group)])  # organism names, removing NA
  
  for (name in namesGroup){
    groupRows <- total.data$lab_group%in%c(name, NA)       # remove cases with organisms not of interest. NAs are uninfected
    uniqueGroupRows <- rmOrganism(total.data[groupRows,])
    survDataByGroup[[name]] <- survData.total[groupRows,][uniqueGroupRows,]          # remove duplicate entries
    
    ## include individual group maximum time
    ## omit for overall max e.g. dataset administrative censoring times
    #survDataByGroup[[name]] <- subdistn(survDataByGroup[[name]])
  }
  
  
  ## Gram positive and negative groups
  gram <- c(-1,1)
  for (i in seq_along(gram)){
    groupRows <- total.data$Gram%in%c(gram[i], NA)       # remove cases with organisms not of interest. NAs are uninfected
    uniqueGroupRows <- rmOrganism(total.data[groupRows,])
    survDataByGroup[[as.character(gram[i])]] <- survData.total[groupRows,][uniqueGroupRows,]          # remove duplicate entries
  }
  
  survDataByGroup  
}
## END FUNCTION ##


splitByGroup <- function(total.data, survData.total, groupCol="lab_group"){
  ## more generic form of splitByOrg()
  ## split data a list of into 2 predefined groups with non-HA-BSI for later HR
  ## other groups removed and do not go in with everything else
  ## call: splitByGroup(total.data, survData.total groupCol="combinedMET")
  
  survDataByGroup <- list()
  groupsCol.inf <- total.data[total.data$infstatus==1, groupCol]
  groupsCol.inf[is.na(groupsCol.inf)] <- "no_test"    # replace NAs with group name
  total.data[total.data$infstatus==1, groupCol] <- groupsCol.inf 
  namesGroup <- as.character(unique(groupsCol.inf))   # the group names
  
  for (name in namesGroup){
    ## remove cases with groups not of interest
    groupRows <- (total.data[,groupCol]==name & total.data$infstatus==1) |  total.data$infstatus==0
    groupRows[is.na(groupRows)] <- FALSE
    uniqueGroupRows <- rmOrganism(total.data[groupRows,])
    survDataByGroup[[name]] <- survData.total[groupRows,][uniqueGroupRows,]          # remove duplicate entries
  }
  
  survDataByGroup
}
## END FUNCTION ##


adminCens <- function(total.data){
#   if administratively censored
#   fill-in censoring time and add flag
#   infected and non-infected datasets handled separately  
#   
#   maximum discharge/death date for each sample separately
#   i.e. administrative censoring times (not necessarily equal)
  mix.disdteMax <- max(total.data$hes_disdte[total.data$infstatus==0], na.rm=TRUE)
  inf.disdteMax <- max(total.data$hes_disdte[total.data$infstatus==1], na.rm=TRUE)
    
  # missing discharge/death time & not left hospital yet, logical vector
  inhosp <- total.data$missingTime==0 &
            total.data$hes_dismethdescription%in%c("Not applicable patient still in hospital", "Not known") 
  
  # not infected
  total.data$hes_disdte[total.data$infstatus==0 & inhosp] <- mix.disdteMax
  
  # infected
  total.data$hes_disdte[total.data$infstatus==1 & inhosp] <- inf.disdteMax
  
  # include administratively censored patients with replaced max time as not missing wrt missingTime
  # but with missing event
  total.data$missingTime <- total.data$missingTime | inhosp
  total.data$missingType[inhosp] <- 0
  
  # recalculate updated LOS
  # hospital stay=0 days is (0,1] so replace with 0.5
  # (redundant if only stay>=2 days anyway)
  total.data$adm_diff_dis <- ifelse(total.data$hes_disdte - total.data$hes_admdte==0,
                                         0.5, total.data$hes_disdte - total.data$hes_admdte)
  
total.data
}
## END FUNCTION ##

catStays <- function(total.data){
  ## concatenate same spell records when the
  ## patient is moved within hospital
  ## i.e. discharge=arrival time
  
  uniqueid <- unique(total.data$hes_ID)  # individual patients IDs
  
  progressbar <- txtProgressBar(min = 0, max = length(uniqueid), style = 3)
  counter = 0 
  
  for (i in uniqueid){
    
    counter = counter+1
    setTxtProgressBar(progressbar, counter)
    
    ind <- which(total.data$hes_ID==i)  # identify rows indices of total.data for a single patient
    
    ## TODO ##
    # include admission and discharge description fields in array...
    
    #if ( total.data$disdestdescription == total.data$admisorcdescription %in%
       #   c("NHS other hospital provider - ward for general patients or the younger phys. disabled (from 1 April 2004)",
        #    "NHS other hospital provider - ward for maternity patients or neonates (from 1 April 2004)",
         #   "NHS other hospital provider - high security psychiatric accommodation (from 1 April 2004)")){
      
      ## replace all episodes within a spell with the earliest admission and latest discharge
  
    newTimes <- spellLength(total.data$hes_admdte[ind], total.data$hes_disdte[ind])
    total.data$hes_admdte[ind] <- newTimes$adm
    total.data$hes_disdte[ind] <- newTimes$disch
    
    ## update missing time flag. may be that changed from known to missing
    total.data$missingTime[ind][is.na(total.data$hes_disdte[ind])] <- 0
  }
  close(progressbar)

total.data
}
## END FUNCTION ##

spellLength <- function(adm, disch){
  ## alternative method to calculate
  ## earliest and latest times
  ## for concatenation
  
  for (i in 1:length(adm)){
    ## which records end on start date
    dupl <- adm[i]==disch
    
    ## assign latest discharge
    disch[dupl] <- disch[i]
    
    if(any(dupl, na.rm=TRUE))
      ## assign earliest start
      {adm[i] <- adm[dupl] <- min(adm[dupl], na.rm=TRUE)}
  }

list(adm=adm, disch=disch)
}
## END FUNCTION ##


subdistn <- function(survData){
  # substitute in later times to change risk sets
  # time representation appropriate for subdistribution hazards
  # 
  # Accounts for administrative censoring
  # note that within a specific organism subset the maximum time
  # may be from a patient infected by a different organism
  
  # default entries
  survData$time2disch <- survData$time2death <- survData$time
  
  ## event times of INFECTED database, DISCHARGE alive (censored at death)
  mtime1 <- max(survData[survData$infstatus==1, "time"], na.rm=TRUE)   # maximum infected database time, Dead=TRUE
                                                                      # not actually the administritive censoring time which is bigger than this
  survData$time2disch[survData$infstatus==1 & survData$event & survData$missingType] <- mtime1+1
  
  ## event times of INFECTED database, DEATH (censored at discharge alive)
  survData$time2death[survData$infstatus==1 & !survData$event & survData$missingType] <- mtime1+1
  
  ## event times of UNINFECTED, DISCHARGE alive (censored at death)
  mtime0 <- max(survData[survData$infstatus==0, "time"], na.rm=TRUE)  # maximum uninfected database time, Dead=TRUE
  
  survData$time2disch[survData$infstatus==0 & survData$event & survData$missingType] <- mtime0+1
                            
  ## event times of UNINFECTED, DEATH (censored at discharge alive)
  survData$time2death[survData$infstatus==0 & !survData$event & survData$missingType] <- mtime0+1

survData
}
## END FUNCTION ##


stratSurvData <- function(survData, survDataCol){
##
## DEPRECATED BY split(base)
## stratify sample into each factor group

	Glevels <- nlevels(survDataCol) # number of factor levels
	Gnames <- levels(survDataCol)   # names of factor levels
	survDataStrat <- list()

	for (i in 1:Glevels){
		survDataStrat[[i]] <- survData[survDataCol==Gnames[i],]
		survDataStrat[[i]] <- survDataStrat[[i]][order(survDataStrat[[i]]$infstatus),]	# reorder by inf/non-inf
	}

	names(survDataStrat) <- Gnames
	
survDataStrat
}
## END FUNCTION ##



