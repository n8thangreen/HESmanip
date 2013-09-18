
# convert NULL entries to NAs
NULLtoNA <- function(x) as.numeric(as.character(x))


getSurvData <- function(total.data){
  ## create a simpler survival analysis format array
  ## basically a subset of the total.data columns
  ## with the dates removed so all time origins are equal
   
  survData <-  within(total.data, {
    time <- as.numeric(adm_diff_dis)  	# time from admission -> death/discharge. recalculate since including administrative censoring times
    event <- hes_dischargeddeadoralive=="Death"	            # event of interest is death, discharge is censoring time
    age <- hes_ageatstartofepisode			        	          # age at admission
    gender <- droplevels(hes_sexdescription)                # remove unknown levels
    spectime <- as.numeric(adm_diff_specdate)		            # specimen time for non-infected is NA
    
    agegr <- cut(age, breaks=c(0,2,5,10,18,40,60,100,10000))  # cut-up in to grouped ages intervals
    # agegr <- cut(age, breaks=c(0,5,18,100,10000))			    # from JPIDS paper (2012). to check stats agree
    
    rm(lab_opieid, hes_ID, lab_Specimendate, hes_admdte, hes_disdte, adm_diff_specdate, dis_diff_specdate, hes_dismethdescription, adm_diff_dis,
       hes_sexdescription, hes_ageatstartofepisode, hes_dischargeddeadoralive, lab_organismname, lab_group, Gram)  
  })
  
  survData$age[survData$age>1000] <- 0	# convert <1 yr old code to 0 year old
  
  # re-order age factors
  survData$agegr <- factor(survData$agegr,
                           c("(100,1e+04]", "(0,2]", "(2,5]" ,"(5,10]","(10,18]","(18,40]","(40,60]","(60,100]"))
  
  ## include subdistribution event times
  survData <- subdistn(survData)
  
survData
}
## END FUNCTION ##



artifCens <- function(total.data, censtime=NA){
##   artificially impose censoring at either
##   earliest censoring time or a user defined time
##   censtime: fraction of time between smallest and largest administrative censoring times [0,1]
    
  ## minimax of discharge times
  mix.disdteMax <- max(total.data$hes_disdte[total.data$infstatus==0], na.rm=TRUE)
  inf.disdteMax <- max(total.data$hes_disdte[total.data$infstatus==1], na.rm=TRUE)
  
  artifCens.time <- min(mix.disdteMax, inf.disdteMax)
  
  if(is.numeric(censtime) & !is.na(censtime)){
    artifCens.time <- artifCens.time + (max(mix.disdteMax, inf.disdteMax)-artifCens.time)*censtime 
  }
  
  artifCens.dis <- total.data$hes_disdte>artifCens.time # which patients are to be censored
  
  ## set newly censored patient to have missing event times, still in hospital and no event type indicator
  total.data$hes_dismethdescription[artifCens.dis] <- "Not applicable patient still in hospital"
  total.data$hes_dischargeddeadoralive[artifCens.dis] <- ""
  total.data$hes_disdte[artifCens.dis] <- artifCens.time
  total.data$missingType[artifCens.dis] <- 0
  
  ## furthermore
  ## if spectime (infection) is after censoring time too then consider as an uninfected patient only
  artifCens.spec <- total.data$lab_Specimendate>artifCens.time  # which patients are to be censored
  total.data$infstatus[artifCens.spec] <- 0
  total.data$adm_diff_specdate[artifCens.spec] <- total.data$dis_diff_specdate[artifCens.spec] <- total.data$lab_Specimendate[artifCens.spec] <- NA
  
  total.data$adm_diff_dis <- with(total.data, hes_disdte-hes_admdte)
  
total.data
}
## END FUNCTION ##


tvCox <- function(survData, covariates, type){
  ##
  ## time-dependent hazard array representation
  ## can specify which time to use: cause specific or subdistributions (Fine & Gray)
  ##
  ## survData: survival data subset of total.data with indicators
  ## covariates: e.g. c("age", "agegr", "cath", "surgical", "gender")
  ## type: type of hazard functions to use
  ##         "" - Cause-specific
  ##  	"death" - subdistribution time-to-death 
  ##		"alive" - subdistribution time-to-discharge alive 
  
  
  if (type=="alive"){
    time <- survData$time2disch	# subdistribution discharge times
  } else if (type=="death"){
    time <- survData$time2death  # subdistribution death times
  } else if (type==""){ 
    time <- survData$time}        # raw times
  
  rows.inf <- survData$infstatus==1
  
  survData.inf <- survData[rows.inf,]
  survData.mix <- survData[!rows.inf,]
  
  ## repeated pairs of row ids for infected cases
  ninf <- rep(1:nrow(survData.inf), each=2)
  epsilon <- 0.5
  
  time.inf <- time[rows.inf]   # subset of event times for infected cases
  time.mix <- time[!rows.inf]
  
  originalid <- 0							# original array single-line format patient row id
  tvdata.inf <- data.frame(id=ninf,
                           tstart=NA,tstop=NA,inf=NA,disch=NA,death=NA,
                           survData.inf[ninf, c("infstatus",covariates)])
  
  numRows.new <- nrow(tvdata.inf)
  
  if (numRows.new > 0){ # only for infected cases
    for (i in seq(1, numRows.new,by=2)){		# fill two rows at a time
      
      originalid <- originalid+1
      spectime <- survData.inf$spectime[originalid]
      time.inf.id <- time.inf[originalid]
      
      # 1) admission -> infection
      tvdata.inf[i,"tstart"]<- 0 		                              # admission start time always zero
      tvdata.inf[i,"tstop"] <- ifelse(spectime==0, epsilon, spectime)	# so start<stop due to same day events
      tvdata.inf[i,"disch"] <- 0		                              # discharge or censored
      tvdata.inf[i,"inf"] 	<- 0	                            	  # not infected yet
      tvdata.inf[i,"death"] <- 0	                              	# death or censored
      
      
      # 2) infection -> death/discharge
      tvdata.inf[i+1,"tstart"] <- tvdata.inf[i,"tstop"]     # start next state at the end time of first state
      
      # stop time is the event time death, discharge or administrative censoring
      # to ensure start<stop add epsilon due to same day events
      tvdata.inf[i+1,"tstop"]  <- ifelse(time.inf.id==spectime,
                                         time.inf.id+epsilon,
                                         time.inf.id)
      
      ## TODO ##
      ## what if spectime>time? table(survData$spectime>survData$time2disch)
      # currently replaces with NA by default by R
      
      # flag if theres a recorded event type and is the event type discharge
      # otherwise censored
      tvdata.inf[i+1,"disch"] <- (1-as.numeric(survData.inf$event[originalid])) & survData.inf$missingType[originalid]
      
      tvdata.inf[i+1,"inf"] 	<- 1
      
      # flag if there a recorded event type and is the event death
      # otherwise censored
      # censored at sink i.e not infection
      ## TODO ## why is event sometimes in a list??
      tvdata.inf[i+1,"death"] <- unlist(survData.inf$event[originalid]) & survData.inf$missingType[originalid]		
    }
  }
  
  # construct an equivalent matrix for the uninfected patients
  # include arbitrary unique IDs (otherwise bootstrapper ignores patient)
  # split event into death and discharge and can still recover censored data
  tvdata.mix <- data.frame(id=originalid+10+seq_along(survData.mix$time),
                           tstart=0,
                           tstop=ifelse(time.mix>0, time.mix, epsilon),
                           inf=0,
                           disch= (1-as.numeric(survData.mix[,"event"])) & survData.mix[,"missingType"],
                           death= unlist(survData.mix[,"event"]) & survData.mix[,"missingType"],
                           survData.mix[, c("infstatus",covariates)])
  
  
  tvdata <- rbind(tvdata.inf, tvdata.mix)
  tvdata <- cbind(tvdata, status=tvdata$disch+(2*tvdata$death)) # single vector for censored (0), discharge (1) or death (2)
  
  tvdata <- tvdata[!is.na(tvdata[,"tstop"]),]  	# remove empty rows
  
  ## rename some columns
  names(tvdata)[which(names(tvdata)=="gender")] <- "sex"
  names(tvdata)[which(names(tvdata)=="infstatus")] <- "originf"
  
tvdata
}
## FUNCTION END ##


etmFormat <- function(survData, type){
##
## transform data in to etm format
## include state information


  tvdata <- tvCox(survData, covariates=NULL, type=type)

  ## ensure all infected cases are at top of array
  tvdata <- tvdata[order(tvdata$originf, decreasing=TRUE),]
  
  ## number of rows
  numCases.inf <- sum(tvdata$originf==1)
  numCases.mix <- sum(tvdata$originf==0)

  id <- tvdata$id

  entry <- tvdata$tstart		
  exit  <- tvdata$tstop		# failure/event times

  ## states
  ## assumes that all infected cases are at top of array but should be from tvCox anyway
  from <- c(rep(0:1,length.out=numCases.inf), rep(0,numCases.mix))		
  to   <- c(rep(1:2,length.out=numCases.inf), rep(2,numCases.mix))				# discharge & death combined indicator
  to2  <- tvdata$disch+to			# ifelse(tvdata$disch==1,3,to)	# same as status in tvdata i.e. censoring

	## censoring flag on final state where appropriate
	
  ## row ids
  rinf <- 1:numCases.inf                # infected cases
	rmix <- (max(rinf)+1):nrow(tvdata)    # uninfected cases
  
  ## infected cases censored row ids
	rcens.inf <- tvdata$inf[rinf]==1 & (tvdata$disch[rinf]+tvdata$death[rinf])==0 
  
  ## non-infectd cases censored row ids
	rcens.mix <- (tvdata$disch[rmix] + tvdata$death[rmix])==0 
  
	rcens <- c(rcens.inf, rcens.mix)
  
  # replace final state
	to[rcens] <- to2[rcens] <- "cens"

  age <- tvdata$age
  agegr <- tvdata$agegr
  sex <- tvdata$sex

  #msm.data <- data.frame(id,entry,exit,from,to,to2,sex,age,agegr)
  msm.data <- data.frame(id,from,to,entry,exit)
  
  msm.data <- msm.data[!is.na(msm.data$exit),]		# remove NA times

msm.data
}
## FUNCTION END ##


getGrams <- function(data, sign){
  ##
  ## produces a vector of given sign Gram recorded in data
  
  gramDict <- unique(data$lab_group[data$lab_GRAMEXT==sign])
  gramDict <- gramDict[gramDict!="other"]   # remove 'other' group
  
  gramDict
}
## END FUNCTION ##


capwords <- function(s, strict = FALSE) {
  ##
  ## Capitalise words
  ## strict:
  
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
{s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
## END FUNCTION ##

