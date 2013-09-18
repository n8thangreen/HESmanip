byGroupHR <- function(survData.list, naive=T, causespec=T, subdist=T, tofile=FALSE){
##
## Batch calculate all the Cox rgressions
## for each group in turn
## can select which type of HR to calculate
## survDataByGroup: list of survData arrays

  startTime <- Sys.time()
  namesGroup <- names(survData.list)
  results <- list()

  progressbar <- txtProgressBar(min = 0, max = length(namesGroup), style = 3)
  counter = 0 
  
  for (i in namesGroup){
    
    counter = counter+1
    setTxtProgressBar(progressbar, counter)

    timeindnaive <- timedeptcausespec <- timedeptsubdistn <- NA
    if (naive==T){timeindnaive=hr.naive(survData.list[[i]])}
    if (causespec==T){timedeptcausespec=hr.tv.causespecific(survData.list[[i]])}
    if (subdist==T){timedeptsubdistn=hr.tv.subdistribution(survData.list[[i]])}
    
    results[[i]] <- list(timeindnaive=timeindnaive,
                         timedeptcausespec=timedeptcausespec,
                         timedeptsubdistn=timedeptsubdistn)
  }
  
    
  close(progressbar)
  cat(Sys.time()-startTime)
  
  if(tofile==TRUE){
    sink(".\\output\\byOrganismHR.txt"); results; sink()}

results
}
## FUNCTION END ##


byGroupLOS <- function(survData.list, type, standerr=FALSE){
##
## Batch calculate all the cumulative hazards,
## transition probabilites and LOS
## for each of the organism groups in turn

  startTime <- Sys.time()
  namesGroup <- names(survData.list)
  
  progressbar <- txtProgressBar(min = 0, max = length(namesGroup), style = 3)
  counter = 0 

  tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), as.character(0:2)))	# admission-infection-death/discharge
  tra[1, 2:3] <- TRUE
  tra[2, 3] <- TRUE

  results <- list()

  for (i in namesGroup){
    
    counter = counter+1
    setTxtProgressBar(progressbar, counter)
    
    if (i==""){stop("no group name")}
  	if (nrow(survData.list[[i]])==0){		# catch empty sets
  		results[[i]] <- list(NA, NA, NA)
  	}else{
  		results[[i]] <- listLOS(survData.list[[i]], tra, type, standerr)}
  }

  close(progressbar)
  cat(Sys.time()-startTime)

results
}
## FUNCTION END ##


listLOS <- function(survData, tra, type, standerr=FALSE, trim=0.99){
  ##
  ## Collection of length of stay functions and associated outputs
  ## hardcoded number of bootstrap samples
  ##
  ## trim: upper percentile to trim data at
  ## standerr: flag whether to calculate standard error (time consuming)
  ## call: out <- listLOS(survDataByGroup.rm[["all"]], tra, type="", standerr=FALSE)
  
  ## trim outliers in total sample
  #tsurvData <- trimSurv(survData, trim)
  tsurvData <- survData   # no trimmming

  
	msm.data <- etmFormat(tsurvData, type)		# transform into etm format
  #cat("State table:"); print(statetable.msm(from,id,data=msm.data))  ## TODO ## check
	mvna.data <- mvna::mvna(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens")
	etm.data <- etm::etm(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", s=0)#, t=90)
	clos.data <- etm::clos(etm.data, aw=TRUE)

	if(standerr){
		se <- sqrt(var(boot.clos(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 10)))
		#se <- sqrt(var(boot.clos2(data=data.table(msm.data), state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 20)))  # data.table faster version
	}else {se<-0}

list(mvna.data=mvna.data, clos.data=clos.data, se=se, etm.data=etm.data)
}
## FUNCTION END ##


my.ephi <- function(clos.data){
  ##
  ## expected excess length of stay calculated using
  ## components of cLOS output
  
  ## match the times for the events and weights
  clos.new <- merge(cbind(clos.data$time, clos.data$phi.case, clos.data$phi.control),
                    cbind(clos.data$w.time, clos.data$weights), by=1)
  
  ## weighted difference in LOS
  return((clos.new[,2] - clos.new[,3])%*%clos.new[,4])
}
## FUNCTION END ##
