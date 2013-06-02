\name{bw.gwr.lcr}
\alias{bw.gwr.lcr}
\title{Bandwidth selection for locally compensated ridge GWR (GWR-LCR)}
\description{
This function finds an optimal bandwidth for \link{gwr.lcr} via a cross-validation approach
}
\usage{
bw.gwr.lcr(formula, data, kernel="gaussian",
        lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
        adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{lambda}{option for a globally-defined (constant) ridge parameter. Default is lambda=0, which gives a basic GWR fit}
  \item{lambda.adjust}{a locally-varying ridge parameter. Default FALSE, refers to basic GWR without a local ridge adjustment (i.e. lambda=0, everywhere); if TRUE, use cn.tresh to set the maximum condition number.  For locations with a condition number (for its local design matrix) above this user-specified threshold, a local ridge parameter is found}
  \item{cn.thresh}{maximum value for condition number, commonly set between 20 and 30}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
  Returns the adaptive or fixed distance bandwidth
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{ridge, bandwidth, GWR}

