highlevelHCAI <- function(survData.list){
  
  
	###############################
	# Cox regression hazard ratio #
	###############################
	## for death, discharge or both as outcomes of interest
  ## by each organism group
	## using naive, time-dependent cause-specific and time-dependent subdistribution approaches

  
	output.HR <- byGroupHR(survData.list)
	#output.HR <- byGroupHR(survData.list, naive=F, subdist=F, causespec=T, tofile=FALSE)
	#output.HR.cens <- byGroupHR(survData.list.cens)

  
  
	##################
	# Length of stay #
	##################
	## including bootstrapping for se is very slow

	type <- ""	# "alive", "death"	# causes-specific/subdistn event times

  ## double check: drop entries containing NAs or negative length of stay
  survData.list <- lapply(survData.list, function(x) x[!is.na(x$time) & x$time>0,])
  
	## previously all patients are included in non-infected dataset
  ## remove patients with hospital discharge before first infection time
  ## to align with the infected cases dataset
  survData.list.rm <- lapply(survData.list, function(x) x[x$time>=min(x$spectime[x$infstatus==1]),])  # min=2 days  

  
  ## by groups
	output.LOS <- byGroupLOS(survData.list.rm, type, standerr=TRUE)
  
list(output.HR=output.HR, output.LOS=output.LOS)
}
