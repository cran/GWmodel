\name{gwr.collin.diagno}
\alias{gwr.collin.diagno}
\title{Local collinearity diagnostic for GWR}
\description{
This function diagnoses the local collinearity for a basic GWR model
}
\usage{
gwr.collin.diagno(formula, data, bw, kernel="gaussian",
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \code{bw.gwr.lcr};
            fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
  \item{corr.mat}{Local correlation matrix}
  \item{vifs.mat}{Local Variance inflation factors (VIFs) matrix}
  \item{vdp.idx}{Condition index of variance-decomposition proportions}
  \item{vdp.pi}{Variance-decomposition proportions}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{collinearity,diagnostic, GWR}

