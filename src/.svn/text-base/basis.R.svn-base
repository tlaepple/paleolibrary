#Testplan:
#sorting in das Laden der Dateien einbauen !
#+ - scale, detrend [], Bereich und Punktwahl, plot, summary, cor
#Allgemeines NA Handling 30% Schwelle (später da für Poster nicht nötig)


addhistory<-function(x,newhist)
{
	newhist<-paste(date(),newhist)
	attr(x,"history")<-c(attr(x,"history"),newhist)
	return(x)
}




#Coordinate conversion from 1D<->2D
c1t2<-function(x,nLon)
{ 
      x<-x-1  
	lat<-x%/%nLon+1
	lon<-x%%nLon+1
	return(list(lat=lat,lon=lon))
}

c2t1<-function(lat,lon,nLon)
{
	return(nLon*(lat-1)+(lon))
}

mergeattr <- function(data,source1,source2,newhistory='')
{
	result<-data
	temp1<-attributes(source1)
	temp2<-attributes(source2)
	attr(result,'lat')<-c(temp1$lat,temp2$lat)
	attr(result,'lon')<-c(temp1$lon,temp2$lon)
	attr(result,'name')<-c(temp1$name,temp2$name)
	attr(result,'history')<-c(temp1$history,paste(date(),newhistory))
	return(result)
	
}


copyattr <- function(data,source,newhistory='',cclass=TRUE)
{
	
	temp<-attributes(source)
	attr(data,'lat')<-temp$lat
	attr(data,'lon')<-temp$lon
	attr(data,'name')<-temp$name
	attr(data,'history')<-c(temp$history,paste(date(),newhistory))
	if (cclass) class(data)<-class(source)
	return(data)
	
}



is_pTs <- function(data) (sum(class(data) == 'pTs')>0)
is_pField <- function(data) (sum(class(data) == 'pField')>0)

summary.pTs <- function(x, ...)
{
 temp<-attributes(x)
 print('Proxy timeseries object')
 print(paste('Names: ',paste(temp$name,collapse=' / ')))
 print('History')
 print(temp$history)
 print("")
 cat("Time range: ",min(time(x))," - ",max(time(x)), "N:",length(time(x)),"\n")
 cat("Data range: ",min(x)," - ",max(x),"\n")
}

summary.pField <- function(x, ...)
{
 temp<-attributes(x)
 print('Proxy field object')
 print(paste('Names: ',paste(temp$name,collapse=' / ')))
 print('History')
 print(temp$history)
 print("")
 cat("Time range: ",min(time(x))," - ",max(time(x)), "N:",length(time(x)),"\n")
 cat("Data range: ",min(x)," - ",max(x),"\n")
 print("spatial extent ")
 cat('lat: ',min(temp$lat)," - ",max(temp$lat),"N:",length(temp$lat),"\n")
 cat('lon: ',min(temp$lon)," - ",max(temp$lon),"N:",length(temp$lon),"\n")


}

#remove the points which are only containing NA's
prcompNA.pField <- function(data,nPc=2,center=TRUE,scale=TRUE, ...)
  {
      temp<-attributes(data)
        class(data)<-"matrix"
       dat<-data[,!is.na(colSums(data))]
       result<-prcomp(dat,center=center,scale=scale)

       tm<-matrix(NA,ncol(data),ncol(result$rotation))
       tm[!is.na(colSums(data)),]<-result$rotation

       pc<-pTs(result$x[,1:nPc],time(data),paste("PC",1:nPc,temp$name),c(temp$history,"prcomp"))
       eof<-pField(tm[,1:nPc],1:nPc,temp$lat,temp$lon,paste("EOF",temp$name),c(temp$history,"prcomp"))
       sdev<-result$sdev[1:nPc]
       sdsum<-sum(result$sdev)
  return(list(pc=pc,eof=eof,sdev=sdev,sdsum<-sdsum))

  
    }

#apply a function on fields containing complete NA sets...
na.apply<-function(x,FUN,... )
  {
    index<-!is.na(colSums(x))
     x[,index]<-FUN(x[,index], ...)
   return(x)
  }

 

getlat <- function(data) return(attr(data,"lat"))
getlon <- function(data) return(attr(data,"lon"))
getname <- function(data) return(attr(data,"name"))
gethistory <- function(data) return(attr(data,"history"))




maxpoint <- function(data)
  {
  pos<-which(data==max(data))
  value<-max(data)
  lat<-getlat(data)
  lon<-getlon(data)
  
  pos2d<-c1t2(pos,length(lon))
  return(list(lat=lat[pos2d$lat],lon=lon[pos2d$lon],value=value))
  }

minpoint <- function(data)
  {
  pos<-which(data==min(data))
  value<-min(data)
  lat<-getlat(data)
  lon<-getlon(data)
  
  pos2d<-c1t2(pos,length(lon))
  return(list(lat=lat[pos2d$lat],lon=lon[pos2d$lon],value=value))
  }

#apply FUN(field->scalar) for each timestep and gives back a timeseries
applyspace<-function(data,FUN)
{
     index<-!is.na(colSums(data)) 
   ts<-apply(data[,index],1,FUN)
   return(pTs(ts,time(data),name=getname(data)))
       }


#apply FUN(field->scalar) for each gridbox and gives back a single field
applytime<-function(data,FUN,newtime=NULL)
{
   if (is.null(newtime)) newtime<-mean(time(data))
   field<-apply(data,2,FUN)
   return(pField(field,newtime,getlat(data),getlon(data),name=getname(data)))
}


#return 2D Fields filled with lats and lons

latlonField <- function(data)
{
  lat<-getlat(data)
  lon<-getlon(data)

  nlat<-length(lat)
  nlon<-length(lon)
  
  lon2d<-rep(lon,nlat)
  lat2d<-rep(lat,each=nlon)

  
  return(list(lat2d=lat2d,lon2d=lon2d))
}



schwerpunkt<-function(data)
{
#nicht allgemien !
      
      t<-latlonField(data)
      t$lon2d[t$lon2d>180]<- t$lon2d[t$lon2d>180]-360
      lat<-weighted.mean(t$lat2d,data)
      lon<-weighted.mean(t$lon2d,data)
      if (lon < 0) lon<-lon+360
      return(list(lat=lat,lon=lon))
}

#combine timeseries
cbind.pTs <- function(..., deparse.level = 1)
  {
    print("cbind.pTs")
    result<-cbind.ts(..., deparse.level=deparse.level)
    args <- list(...)
    lat<-NULL
    lon<-NULL
    name<-NULL
    for (a in args)
      {
        lat<-c(lat,getlat(a))
        lon<-c(lon,getlon(a))
        name<-c(name,getname(a))
      }
    return(pTs(result,time(result),lat,lon,name,"cbind"))
  }

scale_space <- function(data)
{
  data[,]<-scale(as.vector(data))
  return(data)
}


rollmean.pTs <- function(x, k, na.pad = TRUE, align = c("center", "left", "right"), ...)
{
  return(applyData(x,rollmean,k,na.pad, align, ...))
}

applyData<-function(x,fun,... )
  {
    x[]<-fun(as.vector(x),... )
    return(x)  }



