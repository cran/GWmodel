\name{mink.matrixview}
\alias{mink.matrixview}
\title{Visualisation of the results from \code{\link{mink.approach}}}
\description{
This function visualises the AICc/CV results from the \code{\link{mink.approach}}.
}
\usage{
mink.matrixview(diag.df, znm=colnames(diag.df)[4], criterion="AIC")
}

\arguments{
  \item{diag.df}{the first part of a list object returned by \code{\link{mink.approach}}}
  \item{znm}{the name of the forth column in diag.df}
  \item{criterion}{the criterion used for distance metric selection in \code{\link{mink.approach}}}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{Minkovski approach, visualization}

