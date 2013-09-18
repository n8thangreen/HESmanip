

long2wide <- function(data.long){
  #
  # convert array from etm long format to "old" wide format
  # data.wide <- long2wide(data.long)
  # data.wide <- long2wide(dat.icu.pneu)
  
  data.wide <- data.frame(adm.id=sort(unique(data.long$id)))
  data.wide[,c("j.01","j.02", "j.03","j.12","j.13","cens")] <- Inf
  if ("time"%in%names(data.long)){
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==1], "j.01"] <- data.long$time[data.long$from==0 & data.long$to==1]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==2], "j.02"] <- data.long$time[data.long$from==0 & data.long$to==2]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==3], "j.03"] <- data.long$time[data.long$from==0 & data.long$to==3]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==1 & data.long$to==2], "j.12"] <- data.long$time[data.long$from==1 & data.long$to==2]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==1 & data.long$to==3], "j.13"] <- data.long$time[data.long$from==1 & data.long$to==3]
    
    ## censored times
    data.wide[data.wide$adm.id%in%data.long$id[data.long$to=="cens"], "cens"] <- data.long$time[data.long$to=="cens"]  
  }
  if ("exit"%in%names(data.long)){
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==1], "j.01"] <- data.long$exit[data.long$from==0 & data.long$to==1]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==2], "j.02"] <- data.long$exit[data.long$from==0 & data.long$to==2]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==0 & data.long$to==3], "j.03"] <- data.long$exit[data.long$from==0 & data.long$to==3]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==1 & data.long$to==2], "j.12"] <- data.long$exit[data.long$from==1 & data.long$to==2]
    data.wide[data.wide$adm.id%in%data.long$id[data.long$from==1 & data.long$to==3], "j.13"] <- data.long$exit[data.long$from==1 & data.long$to==3]
    
    ## censored times
    for (i in data.long$id[data.long$to=="cens"]){
      if(sum(data.long$id==i)>1){
        data.wide[data.wide$adm.id==i, "cens"] <- data.long$exit[data.long$id==i][2]
      }else{
        data.wide[data.wide$adm.id==i, "cens"] <- data.long$exit[data.long$id==i][1]}
    }
  }
  
 data.wide[data.wide$j.12==Inf & data.wide$j.02==Inf &data.wide$j.03==Inf &data.wide$cens==Inf &data.wide$j.13==Inf &data.wide$j.01!=Inf,"cens"]<- data.wide[data.wide$j.12==Inf & data.wide$j.02==Inf &data.wide$j.03==Inf &data.wide$cens==Inf &data.wide$j.13==Inf &data.wide$j.01!=Inf, "j.01"]

data.wide
}




cLOS <- function(x=NA, my.data)
{
  ## This program is published under the terms of the GNU General
  ## Public License. The terms of GNU General Public License can be
  ## obtained at http://www.gnu.org/copyleft/gpl.html. The program
  ## is free software and comes with absolutely no warranty. It may
  ## be redistributed under the conditions of the GNU
  ## General Public License.
  ## Author: Jan Beyersmann, jan@fdm.uni-freiburg.de
  
  require(survival)
  
  ## need x for bootstrapping, e. g.:
  library(bootstrap)
  
  ## result <- bootstrap(x=1:length(los.data[,1]), nboot=50, theta=function(x){y <- cLOS(x); return(y$cLOS)}, func=var)
  ##if(is.na(x[1])){
  ## war bis 19.9.04: 
  if(is.na(x)){
    x <- 1:length(my.data[,1])
  }
  my.data <- my.data[x,]
  
  ## compute variables cens.0 for admissions censored in the
  ## initial state 0 and cens.1 for admissions censored in state 1
  
  my.data$cens.0 <- my.data$cens
  my.data$cens.0[is.finite(my.data$j.01)] <- Inf
  my.data$cens.1 <- my.data$cens
  my.data$cens.1[is.infinite(my.data$j.01)] <- Inf
  
  ## compute `transition matrix' for every jump time
  jump.times <- sort(unique(c(my.data$j.01, my.data$j.02,
                              my.data$j.03, my.data$j.12, my.data$j.13)))
  jump.times <- jump.times[is.finite(jump.times)]
  jump.matrices <- array(0, c(4,4,length(jump.times)))
  
  for(i in 1:length(jump.times)){
    ## compute number of jumps
    jump.matrices[1,2,i] <- length(my.data$j.01[my.data$j.01==jump.times[i]]) ## jump 0->1
    jump.matrices[1,3,i] <- length(my.data$j.02[my.data$j.02==jump.times[i]]) ## jump 0->2
    jump.matrices[1,4,i] <- length(my.data$j.03[my.data$j.03==jump.times[i]])
    # etc.
    jump.matrices[2,3,i] <-length(my.data$j.12[my.data$j.12==jump.times[i]])
    jump.matrices[2,4,i] <-length(my.data$j.13[my.data$j.13==jump.times[i]])
    
    
    ## divide by number at risk (only necessary, if risk set is not empty)
    risk.0 <- length(my.data[,1]) - length(c(my.data$j.01, my.data$j.02, my.data$j.03, my.data$cens.0)
                                           [c(my.data$j.01, my.data$j.02, my.data$j.03, my.data$cens.0)< jump.times[i]])
    
    risk.1 <- length(my.data$j.01[my.data$j.01 < jump.times[i]])- length(c(my.data$j.12, my.data$j.13,my.data$cens.1)
                                                                         [c(my.data$j.12, my.data$j.13, my.data$cens.1)< jump.times[i]])
    
    if(risk.0 > 0) jump.matrices[1,,i]<- jump.matrices[1,,i]/risk.0
    if(risk.1 > 0) jump.matrices[2,,i]<- jump.matrices[2,,i]/risk.1
    
    ## compute diagonal elements; sum over each row should be 1.
    jump.matrices[,,i] <- jump.matrices[,,i] + diag(1,4,4) - diag(apply(jump.matrices[,,i], 1, sum))
  }
  ## end of compute transition matrices
  
  
  
  ## compute expected LOS given the state at every observed
  ## transition time _except for_ the greatest observed time
  ## (which may be a censoring time)
  ## is there a censoring time greater than the last observed
  ## transition time?
  my.times <- sort(unique(c(jump.times, max(my.data$cens[is.finite(my.data$cens)], jump.times))))
  
  los <- matrix(data=rep(my.times,3), ncol=3, byrow=FALSE, dimnames=list(NULL, c("Time", "Given in state 1","Given in state 0")))
  los[length(my.times)-1,2:3] <- rep(max(my.times), 2)
  
  
  ## added my NGreen ##
  ## median length of stay
  med.los <- matrix(data=rep(my.times,3), ncol=3, byrow=FALSE, dimnames=list(NULL, c("Time", "Given in state 1","Given in state 0")))
  
  
  ## last two rows in los already correct. compute the rest,
  ## starting with the last but two row (which corresponds to the
  ## third greatest time in my.times)
  ## will need to temporarily store Aalen-Johansen estimates
  
  aj <<- array(NA, c(4,4,1))
  aj[,,1] <- diag(1,4,4)
  
  ## will need function that does matrix multiplication running
  ## thru the `slices' of array aj 
  
  my.function <- function(x,y){x%*%aj[,,y]}
  
  for(i in (length(my.times)-2):1){
    ## compute time differences
    
    diffs <- diff(my.times[(i+1):length(my.times)])
    
    ## find `transition matrix' that corresponds to my.times[i + 1]...
    
    my.matrix <- jump.matrices[,,length(jump.times[jump.times <= my.times[i + 1]])]
    
    ## ...and multiply it with Aalen-Johansen estimates of the previous loop
    
    aj <- array(apply(X=diag(1:dim(aj)[3]), 1, my.function, x=my.matrix), c(4,4,dim(aj)[3]))
    
    ## LOS given in state 1 at time my.times[i]
    los[,2][i] <- my.times[i+1] + matrix(diffs, nrow=1) %*% matrix(aj[2,2,],ncol=1)
    ## LOS given in state 0 at time my.times[i]
    los[,3][i] <- my.times[i+1] + matrix(diffs, nrow=1) %*% matrix((aj[1,1,] + aj[1,2,]),ncol=1)
    
    
    ## added by NGreen ##
    ## median modification
    ecdf.inf <- aj[2,3,]+aj[2,4,]
    med.min <- min(which(ecdf.inf >= 0.5))
    med.los[,2][i] <- my.times[i + med.min]
    
    ecdf.adm <- aj[1,3,]+aj[1,4,]
    med.min <- min(which(ecdf.adm >= 0.5))
    med.los[,3][i] <- my.times[i + med.min]
    
    
    
    ## stack identity matrix on top for the next loop
    aj <- array(c(diag(1,4,4), aj), c(4,4, (dim(aj)[3] + 1)))
  }
  
  
  ## compute distribution to weight differences in LOS
  ## need waiting time distribution in initial state 0.
  ## create a survival object and fit it.
  ## event: left state 0.
  
  T0 <- Surv(apply(as.matrix(my.data[, c("j.01","j.02", "j.03","cens")]), 1, min), 1-is.finite(my.data$cens.0))
  T0.fit <- survfit(T0~1)
  
  ## Only need `time' and `surv' for non-censoring events:
  T0.fit$time <- T0.fit$time[T0.fit$n.event!=0]
  T0.fit$surv <- T0.fit$surv[T0.fit$n.event!=0]
  
  ## weight by waiting time distribution in initial state 0
  ## need mass for each event time my.weights[i] corresponds to T0.fit$time[i]
  my.weights <- diff(c(0, 1-T0.fit$surv))
  
  ## compute estimate
  estimate <- matrix((los[,2]-los[,3])[is.element(los[,1], T0.fit$time)], nrow=1) %*% matrix(my.weights, ncol=1)
  
  ## added by NGreen ##
  ## median modification
  diff.e <- matrix((med.los[,2]-med.los[,3])[is.element(med.los[,1], T0.fit$time)], nrow=1)
  nNA <- !is.na(diff.e)
  diff.e <- diff.e[nNA]
  estimate.med <- diff.e %*% matrix(my.weights[nNA], ncol=1)
  
  
  
  ## return results
  my.return <- list(cLOS=estimate, cLOS.med=estimate.med, times=jump.times, e.given.1=c(los[,2]), e.given.0=c(los[,3]), e.given.1.med=c(med.los[,2]), e.given.0.med=c(med.los[,3]), weights=my.weights, matrices=jump.matrices)
  
  return(my.return)
  
  
  ### alternative weighting approach. Currently not used ###
  ## 2nd: waiting time distribution in initial state 0
  ## given cause for leaving is state 1
  
  #pr.cause1 <-matrix(c(1, T0.fit$surv[1:length(T0.fit$surv)-1]),nrow=1) %*%
  #matrix(jump.matrices[1,2,][is.element(jump.times,T0.fit$time)], ncol=1)
  
  ## my.weights.1[i] corresponds to T0.fit$time[i]
  ##my.weights.1 <- diag(diag(c(1, T0.fit$surv[1:length(T0.fit$surv)-1])) %*% diag(jump.matrices[1,2,][is.element(jump.times, T0.fit$time)]))/ pr.cause1
  
  #estimate.1 <- matrix((los[,2]-los[,3])[is.element(los[,1],T0.fit$time)], nrow=1) %*% matrix(my.weights.1, ncol=1)
  
  ## 3rd: waiting time distribution in initial state 0 given
  ## cause for leaving are states 2 or 3
  
  #pr.cause23 <- matrix(c(1,T0.fit$surv[1:length(T0.fit$surv)-1]),nrow=1) %*%
  #matrix((jump.matrices[1,3,] + jump.matrices[1,4,])[is.element(jump.times, T0.fit$time)], ncol=1)
  
  ## my.weights.23[i] corresponds to T0.fit$time[i]
  ##my.weights.23 <- diag(diag(c(1, T0.fit$surv[1:length(T0.fit$surv)-1])) %*% diag((jump.matrices[1,3,] + jump.matrices[1,4,])[is.element(jump.times, T0.fit$time)])) /
  ## pr.cause23
  
  #estimate.23 <- matrix((los[,2]-los[,3])[is.element(los[,1],T0.fit$time)], nrow=1) %*% matrix(my.weights.23, ncol=1)
  
  ## return results
  ##my.return <- list(overall=estimate, given.1=estimate.1, given.no1=estimate.23)
  
  #return(my.return)
}
## end of function