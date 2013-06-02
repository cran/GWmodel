\name{glyph.plot}
\alias{glyph.plot}
\title{Multivariate glyph plots of GWPCA loadings}
\description{
This function provides a multivariate glyph plot of GWPCA loadings at each output location.
}
\usage{
glyph.plot(ld,loc, r1=50, add=FALSE,alpha=1,sep.contrasts=FALSE)
}

\arguments{
  \item{ld}{GWPCA loadings returned by \link{gwpca}}
  \item{loc}{a two-column numeric array for providing evaluation locations of GWPCA calibration}
  \item{r1}{argument for the size of the glyphs, default is 50}
  \item{add}{if TRUE, add the plot to the existing window.}
  \item{alpha}{the level of transparency of glyph from function rgb() and ranges from 0 to max (fully transparent to opaque)}
  \item{sep.contrasts}{allows different types of glyphs and relates to whether absolute loadings are used (TRUE) or not}
}
\references{
Harris P, Brunsdon C, Charlton M (2011)
Geographically weighted principal components analysis.
International Journal of Geographical Information Science 25:1717-1736
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{glyph plot, GWPCA}

