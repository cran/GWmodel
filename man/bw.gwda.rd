\name{bw.gwda}
\alias{bw.gwda}
\alias{wqda.cr}
\alias{wlda.cr}
\title{Bandwidth selection for GW Discriminant Analysis}
\description{
A function for bandwidth selection for GW Discriminant Analysis
}
\usage{
bw.gwda(formula, data, COV.gw = T, prior.gw = T, mean.gw = T,
                 prior = NULL, wqda = T, kernel = "gaussian", adaptive
                 = FALSE, p = 2, theta = 0, longlat = F)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame for training, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{COV.gw}{if true, localised variance-covariance matrix is used for GW discriminant analysis; otherwise, global variance-covariance matrix is used}
  \item{mean.gw}{if true, localised mean is used for GW discriminant analysis; otherwise, global mean is used}
  \item{prior.gw}{if true, localised prior probability is used for GW discriminant analysis; otherwise, fixed prior probability is used}
  \item{prior}{a vector of given prior probability}
  \item{wqda}{if TRUE, weighted quadratic discriminant analysis will be applied; otherwise weighted linear discriminant analysis will be applied}
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
}
\value{
An optimum bandwidth
}
\references{
Brunsdon, C, Fotheringham S,  and Charlton, M (2007),
Geographically Weighted Discriminant Analysis, Geographical Analysis 39: 376-396
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{GWDA}

