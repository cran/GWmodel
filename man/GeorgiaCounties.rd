\name{GeorgiaCounties}
\alias{Gedu.counties}
\docType{data}
\title{Georgia counties data (SpatialPolygonsDataFrame)}
\description{
  The Georgia counties data used for Georgia census data.
}
\usage{data(GeorgiaCounties)}
\details{
Variables are from GWR3 file GData_utm.csv.
}
\examples{
data(GeorgiaCounties)
plot(Gedu.counties)
data(Georgia)
coords <- cbind(Gedu.df$X, Gedu.df$Y)
educ.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
plot(educ.spdf, add=TRUE)

}
\keyword{datasets}
