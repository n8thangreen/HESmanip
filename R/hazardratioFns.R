#####################################
#
# Cox regression aggregated functions
#
# Nathan Green
# 11-2012
#
#####################################


hr.naive <- function(survData){
#
# naive (time-independent) Cox regression
# cause-specific hazards
#
# admission -> death (discharge), with discharge (death) as censored times
# NB can't have Gram and group together since one is coarsening of other so leads to singularity

  ## replace NAs with 0
  survData$spectime <- ifelse(is.na(survData$spectime), 0, survData$spectime)
  
  ## covariates string
  cov <- paste(c(names(survData)[1:17], "age", "gender", "infstatus"), collapse= "+")
  
  ## hazard of discharge alive (censored at death)
  ## include censored patients with dead
  #survData$event[survData$missingType==0] <- TRUE		# censor at missing event type (equivalent to death). set as TRUE since taking negative in regression
  #y <- Surv(time, 1-as.logical(event))
  #fit.disch <- summary(coxph(as.formula(paste(y, " ~ ", cov)), data=survData))
  
  fit.disch <- summary(coxph(Surv(time, 1-as.logical(event)) ~age+gender+Tai+other+head+cath+surgical+prem+cancer+
                               hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+infstatus, data=survData))
                             #  age+gender+spectime+cath+surg+prem+cancer, data=survData))

  ## hazard of death (censored at discharged alive
  ## include censored patients with discharge
  survData$event[survData$missingType==0] <- FALSE		# censor at missing event type (equivalent to discharge alive)
  fit.dead <- summary(coxph(Surv(time, as.logical(event)) ~age+gender+Tai+other+head+cath+surgical+prem+cancer+
                              hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+infstatus, data=survData))
                            #  age+gender+spectime+cath+surg+prem+cancer, data=survData))

  # hazard of death or discharge
  fit.both <- summary(coxph(Surv(time, missingType) ~ age+gender+Tai+other+head+cath+surgical+prem+cancer+
                              hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+infstatus, data=survData))
                             # gender+agegr+spectime+cath+surg+prem+cancer, data=survData))

return(list(fit.disch=fit.disch, fit.dead=fit.dead, fit.both=fit.both))
}
## END FUNCTION ##


hr.tv.causespecific <- function(survData){
#
# cause-specific Cox regression
# time-dependent covariate
#
# dummy variable for infection status 
# Convert to split non-infection and infection periods over 2 lines for each infected individual with tvCox()
#

  ## create time-dependent covariate format survival array
  tvdata <- tvCox(survData, type="")	
  
  
  ## TODO ##
  # replace explicit covariate names in formulae with derived string
  ## covariates
  #cov <- paste(names(tvdata)[??], cluster(id), age, sex, inf, collapse= "+")
  
  ## summary() includes extra output information
  
  ## hazard of discharge alive (censored at death and infection)
  #fit.alive.timeonly <- summary(coxph(Surv(tstop-tstart, as.numeric(disch))~inf+cluster(id), data=tvdata))			    # infection time only
  y.alive <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status==1)
  
  fit.alive.timeonly <- summary(coxph(y.alive ~cluster(id)+inf, data=tvdata))  		    # infection time only

#  fit.alive.full <- summary(coxph(as.formula(paste(y.alive, " ~ ", cov)), data=survData))
  fit.alive.full <- summary(coxph(y.alive ~cluster(id)+age+sex+Tai+other+head+cath+surgical+strata(prem)+cancer+
                                    hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+inf, data=tvdata))		# fully adjusted

  ## hazard of death (censored at discharge and infection)
  y.dead <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status==2)
  
  fit.dead.timeonly <- summary(coxph(y.dead ~cluster(id)+inf, data=tvdata))          # infection time only
  
  #  fit.dead.full <- summary(coxph(as.formula(paste(y.dead, " ~ ", cov)), data=survData))
    fit.dead.full <- summary(coxph(y.dead ~cluster(id)+age+sex+Tai+other+head+cath+surgical+strata(prem)+cancer+
                                   hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+inf, data=tvdata))  	# fully adjusted

  ## hazard of death or discharge
  y.both <- Surv(tvdata$tstart, tvdata$tstop, tvdata$status!=0)
  
  fit.both.timeonly <- summary(coxph(y.both ~cluster(id)+inf, data=tvdata))		# infection time only
  
  #  fit.both.full <- summary(coxph(as.formula(paste(y.both, " ~ ", cov)), data=survData))
  fit.both.full <- summary(coxph(y.both ~cluster(id)+age+sex+Tai+other+head+cath+surgical+strata(prem)+cancer+
                                   hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+inf, data=tvdata))	# fully adjusted

  
return(list(fit.alive.timeonly=fit.alive.timeonly, fit.alive.full=fit.alive.full,
			      fit.dead.timeonly=fit.dead.timeonly, fit.dead.full=fit.dead.full,
				    fit.both.timeonly=fit.both.timeonly, fit.both.full=fit.both.full))
}
## END FUNCTION ##


hr.tv.subdistribution <- function(survData){
#
# Cox regression
# time-dependent subdistribution
#

  ## discharge alive
  tvdata <- tvCox(survData, type="alive")
  fit.alive.timeonly <- summary(coxph(Surv(tstart, tstop, status==1)~cluster(id)+inf, data=tvdata))
  fit.alive.full <- summary(coxph(Surv(tstart, tstop, status==1)~cluster(id)+age+sex+Tai+other+head+cath+surgical+strata(prem)+cancer+
                                    hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+inf, data=tvdata))	

  ## hazard of death
  tvdata <- tvCox(survData, type="death")
  fit.dead.timeonly <- summary(coxph(Surv(tstart, tstop, status==2)~cluster(id)+inf, data=tvdata))
  fit.dead.full <- summary(coxph(Surv(tstart, tstop, status==2)~cluster(id)+age+sex+Tai+other+head+cath+surgical+strata(prem)+cancer+
                                   hes_neocaredescription+hes_admimethdescription+hes_admisorcdescription+inf, data=tvdata))	

return(list(fit.alive.timeonly=fit.alive.timeonly, fit.alive.full=fit.alive.full,
            fit.dead.timeonly=fit.dead.timeonly, fit.dead.full=fit.dead.full))
}
## END FUNCTION ##

