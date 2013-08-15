###################################################
#Bandwidth selection using cross-validation
bw.gwda <- function(formula, data, COV.gw = T, prior.gw = T, mean.gw = T,
                 prior = NULL, wqda = T, kernel = "gaussian", adaptive
                 = FALSE, p = 2, theta = 0, longlat = F)
{
  #data must be given as training data
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
    data <- as(data, "data.frame")
  }
  else
    stop("Given training data must be a Spatial*DataFrame or data.frame object")
  ######## variables from formula
  vars <- all.vars(formula)
  grouping.nm <- vars[1]
  expl.vars <- vars[-1]
  p <- length(expl.vars)
  if (p<2)
     stop("Two or more variables shoule be specfied for analysis")
  #x, y from training data
  res1 <- grouping.xy(data, grouping.nm, expl.vars)
  x<- res1$x
  grouping <- res1$y
  lev <- levels(grouping)
  dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
  dp.n<-nrow(data)
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
  }
  else
  {
    upper<-range(dMat)[2]
    lower<-upper/5000
  }
  bw<-NA
  if(wqda)
     bw <- optimize(wqda.cr,lower=lower,upper=upper,maximum=T,
                     x=x, grouping=grouping, dMat=dMat, COV.gw=COV.gw,
                 mean.gw=mean.gw, prior.gw=prior.gw, prior=prior,
                 kernel = kernel, adaptive =adaptive)
  else
     bw <- optimize(wlda.cr,lower=lower,upper=upper,maximum=T,
                     x=x, grouping=grouping, dMat=dMat, COV.gw=COV.gw,
                 mean.gw=mean.gw, prior.gw=prior.gw, prior=prior,
                 kernel = kernel, adaptive =adaptive)
  bw
}

###Correct ratio
wqda.cr <- function(bw, x, grouping, dMat, COV.gw=T,
                 mean.gw=T, prior.gw=T, prior=NULL,
                 kernel = "gaussian", adaptive = FALSE)
{
  if(adaptive) bw <- round(bw)
  wt <- gw.weight(dMat,bw,kernel,adaptive)
  diag(wt) <- 0
  res.df <- try(wqda(x, grouping, x, wt, COV.gw,
                 mean.gw, prior.gw, prior))
  if(!inherits(res.df, "try-error"))
  {
     n.correct <- length(which(grouping == res.df[,"group.predicted"]))
     correct.ratio <- n.correct/nrow(x)
  }
  else
     correct.ratio <- 0
  if(adaptive)
      cat("Adaptive bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
    else
      cat("Fixed bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
  correct.ratio
}

wlda.cr <- function(bw, x, grouping, dMat, COV.gw=T,
                 mean.gw=T, prior.gw=T, prior=NULL,kernel = "gaussian", adaptive = FALSE)
{
  wt <- gw.weight(dMat,bw,kernel,adaptive)
  diag(wt) <- 0
  res.df <- try(wlda(x, grouping, x, wt, COV.gw,
                 mean.gw, prior.gw, prior))
  if(!inherits(res.df, "try-error"))
  {
     n.correct <- length(which(grouping == res.df[,"group.predicted"]))
     correct.ratio <- n.correct/nrow(x)
  }
  else
     correct.ratio <- 0
  if(adaptive)
      cat("Adaptive bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
    else
      cat("Fixed bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
  correct.ratio
}