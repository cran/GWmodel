\name{gwr.generalised}
\alias{gwr.generalised}
\alias{gwr.binomial}
\alias{gwr.binomial.wt}
\alias{gwr.poisson}
\alias{gwr.poisson.wt}
\alias{gwr.fitted}
\title{Generalised GWR models, including Poisson and Binomial options}
\description{
This function implements generalised GWR
}
\usage{
gwr.generalised(formula, data, regression.points, bw, family ="poisson",
 kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F, dMat, cv=T,tol=1.0e-5, 
 maxiter=20)}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by bw.ggwr();fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{family}{a description of the error distribution and link function to
          be used in the model, which can be specified by \dQuote{poisson} or \dQuote{binomial}}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube:  wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar:   wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{cv}{if TRUE, cross-validation data will be calculated}
  \item{tol}{the threshold that determines the convergence of the IRLS procedure}
  \item{maxiter}{the maximum number of times to try the IRLS procedure}
}
\value{
A list of class \dQuote{ggwrm}:
  \item{GW.arguments}{a \link{list} class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a \link{list} class object including the diagnostic information of the model fitting}
  \item{glm.res}{an object of class inheriting from \dQuote{glm} which inherits from the class \dQuote{lm}, see \link{glm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot.}
  \item{CV}{a data vector consisting of the cross-validation data}
}
\references{
Charlton, M, Fotheringham, S, and Brunsdon, C (2007), GWR3.0, \url{http://gwr.nuim.ie/}.

Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
\dontrun{
DM<-gw.dist(dp.locat=coordinates(londonhp))
bw.f1 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM)
res.poisson<-gwr.generalised(BATH2~FLOORSZ, bw=bw.f1,data=londonhp, dMat=DM)
bw.f2 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM,family ="binomial")
res.binomial<-gwr.generalised(BATH2~FLOORSZ, bw=bw.f2,data=londonhp, dMat=DM,
              family ="binomial")
}
}
\keyword{generalised, GWR}
