\name{gwr.cv}
\alias{gwr.cv}
\title{Cross-validation score for a specified bandwidth for basic GWR}
\description{
This function finds the cross-validation score for a specified bandwidth for basic GWR
}
\usage{
gwr.cv(bw, X, Y, kernel,adaptive, dp.locat, p, theta, longlat,dMat)
}

\arguments{
  \item{bw}{bandwidth used in the weighting function;fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{X}{a numeric matrix of independent samples with extra ones vector  for the 1st column}
  \item{Y}{a one-column vector of dependent samples}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{dp.locat}{a two-column numeric array of observation coordinates}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
  \item{CV.score}{cross-validation score}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{cv, GWR}

