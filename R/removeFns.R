findDuplicates <- function(total.data, nCodesInd){
  # remove duplicate same patient entries on same day of admission with
  # missing event type and observed (dead, alive or administrative censoring)
  # because these have the same admission and discharge data as one another
  # making sure to retain all relevant information
  #
  # nCodesInd: the names of the risk factor/comorbidity column names in total.data
  
  uniqueid <- unique(total.data$hes_ID)  # individual patients IDs
  
  progressbar <- txtProgressBar(min = 0, max = length(uniqueid), style = 3)
  counter = 0 
  
  for (i in uniqueid){
    counter = counter+1
    setTxtProgressBar(progressbar, counter)
    
    ## identify unique admission dates for single patient
    ## because they can be the same or different dates
    
    ind <- which(total.data$hes_ID==i)	# identify rows indices of total.data for a single patient
    
    uniquetime <- unique(total.data$hes_admdte[ind])	
    
    for (j in uniquetime){
      ## number of separate hospital visits for a single patient
      ## identified by distinct admission date
      ## identify row indices of same admission times within a patients subset of rows
      ## i.e. which subset of rows are the same visit
      ## implemented like this because someone else may have the same admission date
      
      subrows <- ind[total.data$hes_admdte[ind]==j]						# corresponding rows indices in total.data for given date
      
      
      if (length(subrows)>1){
        ## if theres more than one record for a single hospital visit
        
        ## merge all risk factor flags within same visit records
        ## so when some records are removed don't lose information
        ##TODO## test this. are column names retained in total.data?
        total.data[subrows, nCodesInd] <- matrix(apply(total.data[subrows, nCodesInd], MARGIN=2, FUN=any),
                                                 nrow=length(subrows), ncol=length(nCodesInd), byrow=TRUE)
        
        ## assumption that the retained record has the same description entries as the deleted records
        ## TODO ##
        ## think of a way of amalgamating these if necessary
        #total.data[subrows, "hes_neocaredescription"]
        #total.data[subrows, "hes_dismethdescription"]
        #total.data[subrows, "hes_distdestdescription"]
        
        ## reorder each patient visit in order of specimen time (infection time)
        ## so that when duplicate are removed to leave the first record
        ## the earliest infection time (so largest infected period) is retained
        ## note places NAs last
        
        if (any(!is.na(total.data$lab_Specimendate[subrows]))){
          total.data[subrows,] <- total.data[subrows,][order(total.data$lab_Specimendate[subrows]),]}
        
        
        if (any(is.na(total.data[subrows, "lab_opieid"])) | anyDuplicated(total.data[subrows, "lab_opieid"])){
          # if there is no lab_id (i.e. assumed complete infection data so not infected)
          # or if the lab_id is a duplicate for that visit
          # assumes only one duplicate of lab_id at the moment
          
          total.data <- rmDuplicates(total.data, subrows)} # remove duplicate rows
      }	
    }
  }
  close(progressbar)

total.data
}
## END FUNCTION ##


rmOrganism <- function(data){
  ## need to handle duplicate entries that only differ by organism type or group e.g same adm and disch
  ## so find duplicate entries except for: lab_id, lab_organismname, lab_group, Gram and specimen dates (i.e. LabBase2 fields)
  ## only retain first entry
  ## data: version of total.data
  ## return: logical vector
  
  ## TODO ##
  # currently keep first record in list which should keep record with first Specimendate (they've been reordered)
  # should this include lab_Specimendate as well? Otherwise assumption that all spec dates per patient visit are equal
  
  ## remove columns not to be matched
  data.drop <- subset(data, select=-c(adm_diff_specdate, dis_diff_specdate, lab_Specimendate,
                                      lab_opieid, lab_organismname, lab_group, Gram,
                                      res, susc, combinedMET, VAN, less1.res, less2.res, less3.res, less4.res, less1.sus, less2.sus, less3.sus, less4.sus))
  
  # or conversely, explicitly specify which fields to match on
  #  data.drop <- subset(data, select=c(hes_ID, hes_admdte))#, hes_disdte, hes_dismethdescription, hes_disdestdescription,
  #                                      adm_diff_dis, hes_sexdescription, hes_ageatstartofepisode, 
  #                                      hes_dischargeddeadoralive, hes_digitprovidercode, hes_neocaredescription,
  #                                      hes_admimethdescription, hes_admisorcdescription, infstatus))
  
  ## unique elements and first occurence of duplicates only
  !duplicated(data.drop, fromLast=FALSE)
}


rmDuplicates <- function(total.data, subrows){
  ##
  ## remove duplicate rows in terms of admission (and discharge) times
  ## for a single patient and single hospital visit
  
  ## logical vector of whether to keep a record
  inc.ind <- hos.ind <- rep(TRUE, nrow(total.data))
  
  ## which of the duplicated rows have missing entries (FALSE)
  ## only over single patient, single visit (same admission date) records
  inc.ind[subrows] <- total.data$hes_dischargeddeadoralive[subrows]%in%c("Death","Live")  # non-missing event type
  hos.ind[subrows] <- total.data$hes_dismethdescription[subrows]%in%c("Not applicable patient still in hospital", "Not known")  # administratively censored
  
  
  if (any(inc.ind[subrows])){
    ## when there's at least one death/discharge
    ## remove all other rows including missing death/discharge rows
    ## and keep only first occurence
    
    total.data <- rmAllButOne(inc.ind, subrows, total.data)
  }
  else if (all(hos.ind[subrows])){
    ## if ALL of the records have patient still in hospital
    ## identify which row of subset is the first entry of patient still in hospital
    ## and keep this entry but remove the others
    
    if (any(!is.na(total.data$hes_disdte[subrows]))){
      ## check that there is at least one non-empty discharge date
      ## identify the maximum discharge date
      
      hos.ind[subrows] <- total.data$hes_disdt[subrows] == max(total.data$hes_disdte[subrows], na.rm=TRUE)
    }

    total.data <- rmAllButOne(hos.ind, subrows, total.data)
  }
  
total.data
}
## END FUNCTION ##


rmAllButOne <- function(ind, subrows, total.data){
  ## identify which row of subset is the first entry of ind=TRUE
  ## and keep this entry but remove the others
  ##
  ## ind: T/F for all total.data
  ## subrow: rows for particular patient stay
  
  ind.min <- min(which(ind[subrows]))
  ind[subrows] <- FALSE
  ind[subrows[ind.min]] <- TRUE
  total.data <- total.data[ind,]  # remove all but first in-hospital rows

total.data
}
## END FUNCTION ##


trimSurv <- function(survData, percentile){
  #
  # trim (upper) extreme values for infected and non-infected separately
  #
  qdata.mix <- quantile(survData[survData$infstatus==0,]$time, c(1-percentile, percentile))
  qdata.inf <- quantile(survData[survData$infstatus==1,]$time, c(1-percentile, percentile))
  
  survData <- survData[
    ## one-sided
    #(survData$time<=qdata.mix[2] & survData$infstatus==0) | (survData$time<=qdata.inf[2] & survData$infstatus==1),]
  
    ## one-sided, non-inf only
    (survData$time<=qdata.mix[2] & survData$infstatus==0)| (survData$infstatus==1),]

    ## two-sided
    #(survData$time<=qdata.mix[2] & survData$time>=qdata.mix[1] & survData$infstatus==0) | (survData$time<=qdata.inf[2] & survData$time<=qdata.inf[1] & survData$infstatus==1),]

survData
}
## END FUNCTION ##
