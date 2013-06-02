\name{plot.mcsims}
\alias{plot.mcsims}
\title{Plot the results from the Monte Carlo (randomisation) test of GWPCA}
\description{
This function plots the results from the functions \link{montecarlo.gwpca.1} and \link{montecarlo.gwpca.2}.
}
\usage{
\method{plot}{mcsims}(x, sname="SD of local eigenvalues from randomisations", \dots)
}

\arguments{
  \item{x}{an object of class \dQuote{mcsims}, returned by the function \link{montecarlo.gwpca.1} or \link{montecarlo.gwpca.2}}
  \item{sname}{the label for the observed value on the plot}
  \item{...}{arguments passed through (unused)}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}

\keyword{Monte Carlo, GWPCA}

