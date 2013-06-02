#The spgwr function gwr returns an object of class gwr which contains the parameter estimate and their standard errors. gwr.t.adjust returns the t values for the parameter
#estimates, unadjusted p values, and p values adjusted using (a) Bonferroni, (b) Benjamini-Hochberg, (c) Benjamini-Yekutieli and
#(d) Fotheringham-Byrne procedures.
#Author: MC
#Edited by BL

gwr.t.adjust <- function(gwm.Obj)
{
  hatmatrix <-gwm.Obj$GW.agruments$hatmatrix
  if(!hatmatrix)
     stop("No p-values to be adjusted")
  gwmx <- as.data.frame(gwm.Obj$SDF)
  #colnames(gmx)
  n <- dim(gwmx)[1]
  m <- dim(gwmx)[2]
  vnames<-all.vars(gwm.Obj$GW.agruments$formula)
  vnames[1]<-"Intercept"
  nv <- length(vnames)
  np <- nv
  ntests <- n * np
  enp <- gwm.Obj$GW.diagnostic$enp
  tvals <- as.matrix(gwmx[, (m-1-nv):(m-2)])
  pvals <- 2 * (1 - pt(abs(tvals), ntests))

  bey_pvals <- p.adjust(pvals, "BY", n = ntests)
  beh_pvals <- p.adjust(pvals, "BH", n = ntests)
  bon_pvals <- p.adjust(pvals, "bonferroni", n = ntests)
  dim(bey_pvals) <- c(n,nv)
  dim(beh_pvals) <- c(n,nv)
  dim(bon_pvals) <- c(n,nv)
  #print(bey_pvals)
  vnames <- colnames(tvals)
  colnames(tvals) <- paste(vnames, "_t", sep = "")
  colnames(pvals) <- paste(vnames, "_p", sep = "")
  colnames(bey_pvals) <- paste(vnames, "_by", sep = "")
  colnames(beh_pvals) <- paste(vnames, "_bh", sep = "")
  colnames(bon_pvals) <- paste(vnames, "_bo", sep = "")
  asf_pvals <- pvals * (1 + enp - (enp/ntests))
  asf_pvals[asf_pvals > 1] <- 1
  colnames(asf_pvals) <- paste(vnames, "_fb", sep = "")
  results <- list(t = tvals, p = pvals, by = bey_pvals, fb = asf_pvals,
      bo = bon_pvals, bh = beh_pvals)
  df.res<-data.frame(tvals, pvals, bey_pvals, asf_pvals,bon_pvals, beh_pvals)
  p4s <- proj4string(gwm.Obj$SDF)
  if(is(gwm.Obj$SDF, "SpatialPolygonsDataFrame"))
     polygons<-polygons(gwm.Obj$SDF)
  else
  {
     locat <- coordinates(gwm.Obj$SDF)
     rownames(locat)<-rownames(df.res)
  }
  if (is(gwm.Obj$SDF, "SpatialPolygonsDataFrame"))
  {
     SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=df.res)
  }
  else
     SDF <- SpatialPointsDataFrame(coords=locat, data=df.res, proj4string=CRS(p4s))
  res<-list(results=results, SDF=SDF)
  res
}