LOSboot.plot <- function(survData, rlstns=20, est="mean"){
  ##
  ## rlstns: number of bootstrap samples
  ## est: mean or median estimator
  ##
  ## call: LOSboot.plot(survDataByGroup$all, rlstns=10)

  
  survData <- survData[!is.na(survData$time) & survData$time>0,]
  survData <- survData[survData$time>=min(survData$spectime[survData$infstatus==1]),]  # min=2 days  
  
  cases.inf <- survData.list[[group]]$infstatus==1
  cases.mix <- survData.list[[group]]$infstatus==0
  
  ninf <- sum(cases.inf)
  nmix <- sum(cases.mix)
  
  tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), as.character(0:2)))  # admission-infection-death/discharge
  tra[1, 2:3] <- TRUE
  tra[2, 3] <- TRUE
  
  cLOS <- NULL
  datalist.inf <- datalist.mix <- list()
  
  
  
  for (k in 1:rlstns){
    
    sample.mix.new <- survData[cases.mix,][sample(nmix, replace=T),]
    sample.inf.new <- survData[cases.inf,][sample(ninf, replace=T),]
    
    survData.new <- rbind(sample.mix.new, sample.inf.new)
    
    ## log (skewed) time
    #survData.new$time <- log(survData.new$time)
    #survData.new$spectime <- log(survData.new$spectime)
    
    msm.data <- etmFormat(survData.new, type="")		# transform into etm format
    
    if (est=="mean"){
      mvna.data <- mvna::mvna(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens")
      etm.data <- etm::etm(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", s=0)#, t=90)
      clos.data <- etm::clos(etm.data, aw=TRUE)
      
      ## means ##
      ## vector of bootstrap sample of length of stays
      cLOS[k] <- clos.data$e.phi
      
      ## record each realisation
      datalist.inf[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.case, infstatus=1)
      datalist.mix[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.control, infstatus=0)
    }
    
    if (est=="median"){
      data.wide <- long2wide(msm.data)
      clos.data <- cLOS(my.data=data.wide)
      names(clos.data) <- c("e.phi", "e.phi.med", "times", "phi.case", "phi.control", "phi.case.med", "phi.control.med",  "weights", "matrices")
      
      ## vector of bootstrap sample of length of stays
      cLOS[k] <- clos.data$e.phi.med
      
      ## record each realisation
      datalist.inf[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.case.med, infstatus=1)
      datalist.mix[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.control.med, infstatus=0)
    }
  }
  
  
  ## merge all sets of realisation together
  data.inf <- Reduce(function(x,y) merge(x,y,by=c("time", "infstatus")), datalist.inf)
  data.mix <- Reduce(function(x,y) merge(x,y,by=c("time", "infstatus")), datalist.mix)
  
  data.inf <- data.inf[order(data.inf$time),]
  data.mix <- data.mix[order(data.mix$time),]
  
  data <- rbind(data.inf, data.mix)
  cols <- 2+(1:rlstns)
  
  maxtime <- 50   # x-axis maximum time of conditional stay
  
  #data.all <- cbind(data[,c(1,2)], LOS=apply(data[,cols],1,mean), LOS.lo=apply(data[,cols],1,min), LOS.hi=apply(data[,cols],1,max))
  data.all <- cbind(data[data$time<maxtime,c(1,2)], LOS=apply(data[data$time<maxtime,cols],1,mean), LOS.lo=apply(data[data$time<maxtime,cols],1,mean)-1.96*apply(data[data$time<maxtime,cols],1,sd), LOS.hi=apply(data[data$time<maxtime,cols],1,mean)+1.96*apply(data[data$time<maxtime,cols],1,sd))
  #data.all <- cbind(data[,c(1,2)], LOS=(apply(data[,cols],1,mean)-data$time), LOS.lo=(apply(data[,cols],1,min)-data$time), LOS.hi=(apply(data[,cols],1,max)-data$time))
  
  
  data.all <- cbind(data.all, Status = as.factor(data.all$infstatus))
  
  ##  envelope plot
  x11()
  ggplot(data=data.all, aes(x=time, y=LOS, ymin=LOS.lo, ymax=LOS.hi, fill=Status, linetype=Status)) + geom_line() + geom_ribbon(alpha=0.4)
  
}



LOSboot.plot.wide <- function(data.wide, rlstns=20, est="mean"){
  ##
  ## same as LOSboot.plot but for input data in wide format
  ## rlstns: number of bootstrap samples
  ## est: mean or median estimator
  ##
  ## call: LOSboot.plot.wide(data.wide, rlstns=10)
  
  
  cases.inf <- data.wide$j.01!=Inf
  cases.mix <- data.wide$j.01==Inf
  
  ninf <- sum(cases.inf)
  nmix <- sum(cases.mix)
  
  tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), as.character(0:2)))  # admission-infection-death/discharge
  tra[1, 2:3] <- TRUE
  tra[2, 3] <- TRUE
  
  cLOS <- NULL
  datalist.inf <- datalist.mix <- list()
  
  clos.data.orig <- cLOS(my.data=data.wide)
  
  
  for (k in 1:rlstns){
    
    sample.mix.new <- data.wide[cases.mix,][sample(nmix, replace=T),]
    sample.inf.new <- data.wide[cases.inf,][sample(ninf, replace=T),]
    
    data.wide.new <- rbind(sample.mix.new, sample.inf.new)
        
    clos.data <- cLOS(my.data=data.wide.new)
    names(clos.data) <- c("e.phi", "e.phi.med", "times", "phi.case", "phi.control", "phi.case.med", "phi.control.med",  "weights", "matrices")
    
    if(est=="median"){
      ## vector of bootstrap sample of length of stays
      cLOS[k] <- clos.data$e.phi.med
      
      ## record each realisation
      datalist.inf[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.case.med, infstatus=1)
      datalist.mix[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.control.med, infstatus=0)
    }
    if(est=="mean"){
      ## vector of bootstrap sample of length of stays
      cLOS[k] <- clos.data$e.phi
      
      ## record each realisation
      datalist.inf[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.case, infstatus=1)
      datalist.mix[[k]] <- cbind(time=clos.data$time, phi=clos.data$phi.control, infstatus=0)
    }
  }
  
  
  ## merge all sets of realisation together
  data.inf <- Reduce(function(x,y) merge(x,y,by=c("time", "infstatus")), datalist.inf)
  data.mix <- Reduce(function(x,y) merge(x,y,by=c("time", "infstatus")), datalist.mix)
  
  data.inf <- data.inf[order(data.inf$time),]
  data.mix <- data.mix[order(data.mix$time),]
  
  data <- rbind(data.inf, data.mix)
  cols <- 2+(1:rlstns)
  
  maxtime <- 500   # x-axis maximum time of conditional stay
  
  #data.all <- cbind(data[,c(1,2)], LOS=apply(data[,cols],1,mean), LOS.lo=apply(data[,cols],1,min), LOS.hi=apply(data[,cols],1,max))
  data.all <- cbind(data[data$time<maxtime,c(1,2)], LOS=apply(data[data$time<maxtime,cols],1,mean), LOS.lo=apply(data[data$time<maxtime,cols],1,mean)-1.96*apply(data[data$time<maxtime,cols],1,sd), LOS.hi=apply(data[data$time<maxtime,cols],1,mean)+1.96*apply(data[data$time<maxtime,cols],1,sd))
  #data.all <- cbind(data[,c(1,2)], LOS=(apply(data[,cols],1,mean)-data$time), LOS.lo=(apply(data[,cols],1,min)-data$time), LOS.hi=(apply(data[,cols],1,max)-data$time))
  ## TODO ## data.all <- cbind(data[data$time<maxtime,c(1,2)], LOS=clos.data.orig, LOS.lo=clos.data.orig-1.96*apply(data[data$time<maxtime,cols],1,sd), LOS.hi=clos.data.orig+1.96*apply(data[data$time<maxtime,cols],1,sd))
  
  
  data.all <- cbind(data.all, Status = as.factor(data.all$infstatus))
  
  ##  envelope plot
  x11()
  ggplot(data=data.all, aes(x=time, y=LOS, ymin=LOS.lo, ymax=LOS.hi, fill=Status, linetype=Status)) + geom_line() + geom_ribbon(alpha=0.4)+ ylim(0, 75)#+ylim(0,40)
  
}
