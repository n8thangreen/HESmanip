getHES <- function(dir){
  #
  # reads-in each data file in directory
  # and combines by row into a single array
  # dir: directory containing only STATA .dta files to be read-in
  
  combinedData <- NULL
  filenames <- list.files(dir)  # gets all file names in directory
  
  indivData <- lapply(1:length(filenames), function(x) read.dta(paste(dir, "/", filenames[x], sep="")))
  combinedData <- do.call(rbind, indivData) # combined ages HES dataset
  
  combinedData
}
## END FUNCTION ##