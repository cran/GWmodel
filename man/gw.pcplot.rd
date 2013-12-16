\name{gw.pcplot}
\alias{gw.pcplot}
\title{Geographically weighted parallel coordinate plot for investigating multivariate data sets}
\description{
This function provides a geographically weighted parallel coordinate plot for investigating a multivariate data set.  It has an option that weights the lines of the plot with increasing levels of transparency, according to their observation's distance from a specified focal/observation point.  This plot can be used to identify outliers.
}
\usage{
gw.pcplot(data,vars,focus,bw,ylim=NULL,ylab="",fixtrans=FALSE, p=2, theta=0, 
          longlat=F,dMat,...) 
}

\arguments{
  \item{data}{ a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{vars}{a vector of variable names to be evaluated}
  \item{focus}{an integer, indexing to the observation point}
  \item{bw}{bandwidth used in the weighting function;fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{ylim}{the y limits of the plot}
  \item{ylab}{a label for the y axis}
  \item{fixtrans}{if TRUE, the transparency of the neighbouring observation plot lines increases with distance; If FALSE a standard (non-spatial) parallel coordinate plot is returned.}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{...}{other graphical parameters, (see \link{par})}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\references{
Harris P, Brunsdon C, Charlton M, Juggins S, Clarke A (2013) 
Multivariate spatial outlier detection using robust geographically weighted methods.  
Mathematical Geosciences DOI 10.1007/s11004-013-9491-0

}
\keyword{GW, PCP}

