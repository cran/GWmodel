\name{writeGWR}
\alias{writeGWR}
\title{Write the GWR results}
\description{
This function writes the calibration result of function \link{gwr.basic} to a text file and shape files}
\usage{
writeGWR(x,fn="GWRresults")}
\arguments{
  \item{x}{an object of class \dQuote{gwrm}, returned by the function \link{gwr.basic}}
  \item{fn}{file name for the written results, by default the output files can be found in the working directory, \dQuote{GWRresults.txt}, \dQuote{GWRresults(.shp, .shx, .dbf)}}
}
\seealso{\code{\link{writeGWR.shp}}}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\note{
  The projection file is missing for the writen shapefiles.
}
\keyword{gwr, write}