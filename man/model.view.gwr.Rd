\name{model.view.gwr}
\alias{model.view.gwr}
\title{Visualise the GWR models from \code{\link{model.selection.gwr}}}
\description{
This function visualises the GWR models from \code{\link{model.selection.gwr}}.
}
\usage{
model.view.gwr(DeVar, InDeVars, model.list)
}
\arguments{
  \item{DeVar}{dependent variable}
  \item{InDeVars}{a vector of independent variables for model selection}
  \item{model.list}{a list of all GWR model tried in \code{\link{model.selection.gwr}}}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}

\seealso{\code{\link{model.selection.gwr}}, \code{\link{model.sort.gwr}}}
\examples{
data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
DeVar<-"PURCHASE"
InDeVars<-c("FLOORSZ","GARAGE1","BLDPWW1","BLDPOSTW")
model.sel<-model.selection.gwr(DeVar,InDeVars, data=londonhp,dMat=DM,bw=5000)
model.list<-model.sel[[1]]
model.view.gwr(DeVar, InDeVars, model.list=model.list)
}
\keyword{model, view}

