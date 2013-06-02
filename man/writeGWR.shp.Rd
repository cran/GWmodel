\name{writeGWR.shp}
\alias{writeGWR.shp}
\title{Write GWR results as shape files}
\description{
This function writes a spatial data frame of the calibration result of function \link{gwr.basic} as shape files}
\usage{
writeGWR.shp(x,fn="GWRresults")
}
\arguments{
  \item{x}{an object of class \dQuote{gwrm}, returned by the function \link{gwr.basic}}
  \item{fn}{file name for the written results, by default the output files can be found
            in the working directory, \dQuote{GWRresults(.shp, .shx, .dbf)}}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\seealso{\code{\link{writeGWR}}}
\note{
  The projection file is missing for the written shapefiles.
}
\keyword{gwr, shapefile}
