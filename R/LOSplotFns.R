byGroupHR.boxplot <- function(res2, model, orgs=TRUE){
  ##
  ## orgs: group by organism type or not
  ## call: byorganismHR.boxplot(res2, "dischtime")
  
  ## Report-friendly organism names
  ## full list or subset of most frequent
  ## make check order matches
  
  if(orgs==TRUE){
    ## italicise parts of names
    namesGroup.short <- c("CoNS",
                          "Other",
                          "Non-p. Streptococci",
                          expression(italic("S. aureus")),
                          expression(paste(italic("Enterococcus"), plain(spp.))),            
                          expression(italic("E. Coli")),
                          expression(paste(italic("Klebsiella"), plain(spp.))),
                          expression(paste(italic("Enterobacter"), plain(spp.))),
                          "Gram-negative", "Gram-positive", "All")
    
    #namesGroup.short <- capwords(res2$organism, strict = TRUE)
    
    neworder <- c(11,10,9,1,6,8,5,7,3,4,2)
  }else{
    namesGroup.short <- res2$group
    neworder <- 1:nrow(res2)}
  
  
  ## find column corresponding to model mean
  ## first match
  col <- grep(model, names(res2))
  
  ## plot title
  if (model=="atime HR"){
    model <- "Discharge with time to infection only"
  }
  if (model=="afull HR"){
    model <- "Discharge with full set of covariates"
  }
  if (model=="dtime HR"){
    model <- "Death with time to infection only"
  }
  if (model=="dfull HR"){
    model <- "Death with full set of covariates"
  }
  
  # or alternatively remove title
  # model <- ""
  
  
  wmeans <- as.numeric(res2[,col])[neworder]  # all HR
  liw <- as.numeric(res2[,col+1])[neworder]   # all lower 95% CI
  uiw <- as.numeric(res2[,col+2])[neworder]   # all upper 95% CI
  
  par(mar=c(10, 4, 2, 0.5))	# increase the bottom margin for organism names
  plotCI(barplot(wmeans, col="blue", ylim=c(0,min(10,max(uiw, rm.na=TRUE))), names=namesGroup.short[neworder], las=2, main=model, ylab="Hazard ratio"), wmeans,
         uiw=uiw-wmeans, liw=wmeans-liw, add=TRUE)
  abline(h=1, lty=2, col="red")		# hazard ratio=1
}
## END OF FUNCTION ##



HRboxplot.batch <- function(res2){
##
## disharge time-only/full & death time-only/full 
## individual box plots in new windows
## write to pdf file
## call: HRboxplot.batch(res2)
  
  ## which columns of different event type and covariate subsets
  HRcol <- grep("HR", names(res2))
  HRnames <- names(res2)[HRcol]
  
	for (i in HRnames){
		#x11()
		#pdf(paste(".\\output\\",i,"boxplots.pdf", sep=""))
		byGroupHR.boxplot(res2, i)
    #dev.off()
	}
}
## END OF FUNCTION ##



byGroupLOS.boxplot <- function(output.LOS, offset=0, ADD=FALSE, axis="s", col="black", orgs=TRUE){
  ##
  ## box plots for each group with error bars
  ## offset: shift plot along axis
  ## ADD: overlay multiple plots
  ## orgs: group by organism type or not
  
  wmeans <- se <- NA
  numGroup <- length(output.LOS)
  
  if (orgs==TRUE){
    ## for organism groupings
    ## double-check order of names is correct
    ## italicise parts of names
    namesGroup.short <- c("All",
                          "Gram-positive",
                          "Gram-negative",
                          "CoNS",
                          expression(paste(italic("Enterococcus"), plain(spp.))),   
                          expression(italic("S. aureus")),
                          "Other (Gram-positive)",
                          "Other (Gram-negative)",
                          expression(italic("E. Coli")),
                          "Non-p. Streptococci",
                          expression(paste(italic("Klebsiella"), plain(spp.))),
                          expression(paste(italic("Enterobacter"), plain(spp.))))
    neworder <- c(1,2,3,4,9,12,5,11,10,6,7,8)
  }else{
    ## for other groupings
    namesGroup.short <- names(output.LOS)
    neworder <- 1:numGroup
  }  
  
  ## collect individual group means and standard deviations
  for (i in 1:numGroup){
    if (is.null(output.LOS[[i]]$clos.data$e.phi)){wmeans[i]<-se[i]<-0}
    else{wmeans[i] <- output.LOS[[i]]$clos.data$e.phi
         se[i] <- output.LOS[[i]]$se}
  }
  
  
  wmeans <- unlist(wmeans)[neworder]
  se <- se[neworder]
  #se <- 0.5		# dummy test value
  
  
  
  par(mar=c(5, 9, 2, 0.5))  # increase the left margin for organism names
  bp <- boxplot(0, plot = FALSE, color="yellow") 
  bp$stats <- rbind(wmeans-(2.33*se), wmeans-(1.96*se), wmeans, wmeans+(1.96*se), wmeans + (2.33*se))
  rownames(bp$stats)<-c("","","","","")
  bp$n <- rep(20, numGroup)
  bp$conf <- rbind(wmeans - (2.33*se), wmeans + (2.33*se))
  bp$names <- namesGroup.short[neworder]
  bxp(bp, horizontal=T, xlab="Excess length of stay (days)", boxwex = 0.25, las=1, add=ADD, at=(1:numGroup)+offset, border=col, yaxt=axis)
  abline(v=0, col="red", lty=2)
  
  
  ## Alternative boxplot formats ##
  #forestplot(labeltext=cbind(namesGroup.short), nn=0, zero=0, mean=wmeans, lower=wmeans-se, upper=wmeans+se,
  #  xlab="Time (days)")#, align=NULL) # box size 1/se^2
  
  #plotCI(wmeans,uiw=5, ylab="")		# simple error bar plots
  #abline(h=0, lty=2)
  
  #par(mar=c(10, 4, 2, 0.5))  # increase the bottom margin for organism names
  #plotCI(
  # barplot(wmeans, col="gray", names=namesGroup.short, las=2, main="LOS by organism group with se")#, horiz=TRUE)#, ylim=c(-2,min(20,max(wmeans+se, na.rm=TRUE))))#,
  #  wmeans, uiw=se, add=TRUE)
  
}
## END OF FUNCTION ##


multistateplots.org <- function(output.LOS){
  ##
  ## produce expected length of stay plots
  ## output to file
  
  for (i in names(output.LOS)){
    #pdf(paste(".\\output\\",i,"expLOS.pdf", sep=""))
    #x11()
    plot(output.LOS[[i]]$clos.data, xlim=c(0,40), ylim.e=c(0,90), main=i)
    #dev.off()
    #x11()
    #pdf(paste(".\\output\\",i,"mvnaplot.pdf", sep=""))
    plot(output.LOS[[i]]$mvna.data, xlim=c(0,90), main=i)
    #dev.off()
  }
}
## END FUNCTION ##


trim.boxplot <- function(survData, trim){
  ##
  ## LOS by trimmed mean percentile
  ## trim: vector of percentile points

  for (i in 1:length(trim)){
    out.trim[[i]] <- listLOS(survData, tra, type, standerr=TRUE, trim=trim[i])
    names(out.trim) <- trim
  }
  byGroupLOS.boxplot(out.trim)
}
## END FUNCTION ##


plot.FollowUpChart <- function(id,start,stop,col=1,lty=1,lwd=2,xlim,xlab="",ylab="",...){
  ## plot regpresentation of in-hospital patient journey
  ## admission -> (infection ->) death/discharge
  ## plot.FollowUpChart(data$id, start=data$tstart, stop=data$tstop,
  ##                    xlab="Time since admission (days)", ylab="HA-BSI patients", xlim=c(0,250))
  
  if (missing(xlim)) xlim <- c(min(start),max(stop))
  N <- length(start)
  plot(0,0,type="n",ylim=c(1,N/2),xlim=xlim,xlab="",ylab="",...)  #xlab=xlab,ylab=ylab,xaxt="n",...)
  nix <- lapply(1:N,function(i){
    segments(0,id[i],start[i],id[i] ,col='lightgray',lty=lty,lwd=lwd)  
    segments(start[i],id[i],stop[i],id[i] ,col=col,lty=lty,lwd=lwd)
  })
}
## END FUNCTION ##


plot.FollowUpChartBatch <- function(survData.list, orgs=TRUE){
  ## Creates a single figure of multiple plots by group
  ##
  ## orgs: group by organism type or not
  ## call: plot.FollowUpChartBatch(survDataByGroup)
  
  if(orgs==TRUE){
    ## for organism groupings
    ## double-check ordering is correct
    namesGroup <- c("All",
                    "Gram-positive",
                    "Gram-negative",
                    "CoNS",
                    expression(paste(italic("Enterococcus"), plain(spp.))),   
                    expression(italic("S. aureus")),
                    "Other (Gram-positive)",
                    "Other (Gram-negative)",
                    expression(italic("E. Coli")),
                    "Non-p. Streptococci",
                    expression(paste(italic("Klebsiella"), plain(spp.))),
                    expression(paste(italic("Enterobacter"), plain(spp.))))
    ## subset & reorder groups
    numgroups <- c(1,4,5,6,9,10,11,12)
  }else{  
    ## for other groupings
    namesGroup <- names(survData.list)
    numgroups <- 1:length(namesGroup)
  }  
  
  par(mar=c(5, 6, 2, 1))
  par(oma = c(3, 3, 3, 0))
  par(mai = c(0.5, 0.6, 0.3,0.3))
  
  if(length(numgroups)<6){par(mfrow=c(length(numgroups),1))}else{par(mfrow=c(4,2))}
  
  ## same non-HA-BSI median for all i
  mix.med <- median(survData.list[[1]]$time[survData.list[[1]]$infstatus==0 & survData.list[[1]]$time>=2])
  ## K-M median better but need to extract from survfit object
  #   surv.mix <- Surv(survData.list[[1]]$time[survData.list[[1]]$infstatus==0 & survData.list[[1]]$time>=2],
  #                    survData.list[[1]]$missingType[survData.list[[1]]$infstatus==0 & survData.list[[1]]$time>=2]==1)
  #   fit <- survfit(surv.mix~1)
  #   
  for (i in numgroups){
    ## TODO ##
    ## calculate median times before split by infection time
    inf.med <- median(survData.list[[i]]$time[survData.list[[i]]$infstatus==1])
    spec.med <- median(survData.list[[i]]$spectime, na.rm=TRUE)
    
    data <- tvCox(survData.list[[i]][order(survData.list[[i]]$time),], type="")
    
    ## as proportion of total LoS
    #survData.list[[i]]$spectime <- survData.list[[i]]$spectime/survData.list[[i]]$time
    #survData.list[[i]]$time <- 1
    #data <- tvCox(survData.list[[i]][order(survData.list[[i]]$spectime),], type="")
    
    data <- data[data$originf==1,]
    par(cex.main=2)
    plot.FollowUpChart(data$id, start=data$tstart, stop=data$tstop,
                       xlab="Time since admission (days)", ylab="HA-BSI patients",
                       main=namesGroup[i], xlim=c(0,250))
    points(data$tstop,data$id, pch='l')
    points(data$tstop[data$disch==1], data$id[data$disch==1], pch=19, cex=1)
    points(data$tstop[data$death==1], data$id[data$death==1], pch=8, cex=1, lwd=2, col="red")
    
    
    abline(v=0, col="lightgray")
    abline(v=mix.med)
    abline(v=inf.med, lty=2, lwd=2)
    abline(v=spec.med, lty=3, lwd=2)
  } 
  
  mtext("Time since admission (days)", side=1, outer=T, cex=1.5)
  mtext("HA-BSI patients", side=2, outer=T, cex=1.5)
  #title(main=deparse(substitute(survData.list)), outer=TRUE)
}
## END FUNCTION ##