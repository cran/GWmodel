###Write GWR results
##Author: Binbin Lu

gwr.write<-function(x,fn="GWRresults")
{
   if(!inherits(x, "gwrm")) stop("It's not a gwm object")
   fn1<-paste(fn,".txt",sep="")
   #fn2<-paste(fn,".csv",sep="")
   sink(fn1)
   print.gwrm(x)
   sink()
  # writeGWR.csv(x, fn=fn2)
   gwr.write.shp(x,fn=fn)
   invisible(x)
}
gwr.write.shp<-function(x,fn="GWRresults")
{
   if(!inherits(x, "gwrm")) stop("It's not a gwm object")
   SDF<-x$SDF
   st_write(st_as_sf(SDF), paste(fn, ".shp",sep=""),driver = "ESRI Shapefile",delete_dsn =T)
   invisible(SDF)
}

# This version of this function is kept to make the code work with the early versions of GWmodel (before 2.0-1)
writeGWR<-function(x,fn="GWRresults")
{
   if(!inherits(x, "gwrm")) stop("It's not a gwm object")
   fn1<-paste(fn,".txt",sep="")
   #fn2<-paste(fn,".csv",sep="")
   sink(fn1)
   print.gwrm(x)
   sink()
  # writeGWR.csv(x, fn=fn2)
   gwr.write.shp(x,fn=fn)
   invisible(x)
}

# This version of this function is kept to make the code work with the early versions of GWmodel (before 2.0-1)
writeGWR.shp<-function(x,fn="GWRresults")
{
   if(!inherits(x, "gwrm")) stop("It's not a gwm object")
   SDF<-x$SDF
   st_write(st_as_sf(SDF), paste(fn, ".shp",sep=""),driver = "ESRI Shapefile",delete_dsn =T)
   invisible(SDF)
}