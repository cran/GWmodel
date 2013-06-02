#####Mont Carlo simulation
##Change the positions of a sequence randomly
#Author: Binbin Lu

montecarlo.gwr<-function(formula, data = list(),dMat=NULL,nsims=99, kernel="gaussian",adaptive=F, bw)
{
  ##Extract the model data frame
  this.call <- match.call()
  if (!is.numeric(dMat)||!is(dMat, "matrix"))
      stop("Distance matrix(dMat) has to be specified correctly")
  else if (is.null(dMat))
      dMat <- gw.dist(dp.locat=coordinates(data))  
  if (!is.null(data))
  {
    if (is(data, "Spatial"))
       data <- as(data, "data.frame")
    else
    {
      if (!is(data, "data.frame"))
         stop("Given regression data must be data.frame or Spatial*DataFrame")
    }
  }
  else stop("No regression data frame is avaiable!")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)

    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
  var.n<-ncol(x)
  dp.n<-dim(dMat)[1]
  if (missing(bw))
      stop("Bandwidth must be given for non-adaptive weights")
  if (adaptive)
  {
    stopifnot(is.numeric(bw))
    stopifnot((bw >= 0))
    stopifnot((bw <= 1))
  }
  else
  {
    stopifnot(is.numeric(bw))
    stopifnot((bw > min(dMat)))
  }
  bandwidth<-bw
 #####Calibrate the original GWR model
  betas <- matrix(nrow=dp.n, ncol=var.n)
  Var_betas<-matrix(nrow=nsims+1, ncol=var.n)
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw.reg(x,y,W.i,hatmatrix=F,i)
    betas[i,]<-gw.resi[[1]]
  }
  for (j in 1:var.n)
  {
    Var_betas[1,j]<-var(betas[,j])
  }
 ###################### Random sequence
  sq<-1:length(y)
  for(k in 1:nsims)
  {
    #randsq<-randomSQ(sq)
    #randx<-x[randsq,]
    #randy<-y[randsq]
    mcs <- sample(dp.n)
		dMat[mcs,]<-dMat[1:dp.n,]
		dMat[,mcs]<-dMat[,1:dp.n]
    betas <- matrix(nrow=dp.n, ncol=var.n)
    for (i in 1:dp.n)
    {
      dist.vi<-dMat[,i]
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      gw.resi<-gw.reg(x,y,W.i,hatmatrix=F,i)
      betas[i,]<-gw.resi[[1]]
    }
    for (j in 1:var.n)
    {
      Var_betas[k+1,j]<-var(betas[,j])
    }
  }
  ######Compute the p-values
  p.values<-numeric(var.n)
  for (j in 1:var.n)
  {
     var_betaj<-Var_betas[,j]
     indx<-sort(var_betaj, index.return=T)$ix
     p.values[j]=1-(which(indx==1))/(nsims+1)
  }
  pmat<-matrix(p.values, ncol=1)
  colnames(pmat)<-c("p-value")
  rownames(pmat)<-colnames(x)
  cat("\nTests based on the Monte Carlo significance test\n\n")
  printCoefmat(pmat)
	invisible(pmat)
}
