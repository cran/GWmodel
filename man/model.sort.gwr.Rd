\name{model.sort.gwr}
\alias{model.sort.gwr}
\title{Sort the results of the GWR model selection function \code{\link{model.selection.gwr}}.}
\description{
Sort the results of the GWR model seclection function \code{\link{model.selection.gwr}}
}
\usage{
model.sort.gwr(Sorting.list , numVars, ruler.vector)
}
\arguments{
  \item{Sorting.list}{a list returned by function \code{\link{model.selection.gwr}}}
  \item{numVars}{the number of independent variables involved in model selection}
  \item{ruler.vector}{a numeric vector as the sorting basis}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\note{
  The function sorts the results of model selection within invidual levels.
}
\seealso{\code{\link{model.selection.gwr}}, \code{\link{model.view.gwr}}}
\keyword{model, sort}

