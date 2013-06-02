\name{gwr.t.adjust}
\alias{gwr.t.adjust}
\title{Adjust p-values for multiple hypothesis tests in basic GWR}
\description{
Given a set of p-values from the t-tests in GWR outputs, returns p-values adjusted using one of several methods, including (a) Bonferroni, 
(b) Benjamini-Hochberg, (c) Benjamini-Yekutieli and (d) Fotheringham-Byrne procedures..
}
\usage{
gwr.t.adjust(gwm.Obj)
}
\arguments{
  \item{gwm.Obj}{an object of class \dQuote{gwrm}, returned by the function \link{gwr.basic}}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289C300. 

Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165C1188. 

Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65C70. 

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383C386. 

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800C803. 

Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561C576. (An excellent review of the area.) 

Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494C504. 

Sarkar, S., and Chang, C. K. (1997). Simes' method for multiple hypothesis testing with positively dependent test statistics. Journal of the American Statistical Association, 92, 1601C1608. 

Wright, S. P. (1992). Adjusted P-values for simultaneous inference. Biometrics, 48, 1005C1013. (Explains the adjusted P-value approach.) 

Byrne, G., Charlton, M. and Fotheringham, S., 2009. Multiple dependent hypothesis tests in geographically weighted regression. In: Lees, B. and Laffan, S. eds. 10th International conference on geocomputation. Sydney.
}
\keyword{gwr, p-values, adjustment}
