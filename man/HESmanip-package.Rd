\name{HESmanip-package}
\alias{HESmanip-package}
\alias{HESmanip}
\docType{package}
\title{Hospital length of stay data preprocessing and analysis}
\description{Hospital length of stay data preprocessing and analysis}
\details{
\tabular{ll}{
Package: \tab HESmanip\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-06-19\cr
License: \tab What license is it under?\cr
}
This package was originally devised and used to analyse a linked LabBase 2 and  Hospital Episode Statistics data set for paediatric patients in England.
Many of the field names and labelling in general are specific to this dataset. However, the fundamental preprocessing and batch calculation approached can be modified and applied to more general data sets.
}
\author{
N Green

Maintainer: N Green <nathan.green@phe.gov.uk>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
#Rprof("highlevelHCAI.out")   # profile function for testing
output <- highlevelHCAI(survDataByGroup)
#Rprof(NULL)
}
