
latmean<-function(data,na.rm=FALSE)
{
	temp<-attributes(data)
	stime<-time(data)
	nTime<-length(time(data))
	nLat<-length(temp$lat)
  	nLon<-length(temp$lon)
	data<-as.matrix(data)
	dim(data)<-c(nTime,nLon,nLat)
	zmeans<-apply(data,c(1,3),mean,na.rm=na.rm)
	return(zmeans)
}

hovmoeller<-function(data,refperiod=c(start(data)[1],end(data)[1]),xlab="time",ylab="latitude",FUN=contour(time(data),getlat(data),zmeans.anomaly,add=T), ...)
{
	zmeans<-latmean(data)
	refmeans<-latmean(window(data,refperiod[1],refperiod[2]))

	zmeans.anomaly<-zmeans-rep(colMeans(refmeans),each=dim(zmeans)[1])
      ##Reference period
  	 palette=colorRampPalette(c("violetred4","blue","steelblue", "lightgreen","white", "yellow","orange","red","brown"))
                         	filled.contour.own(time(data),getlat(data),zmeans.anomaly,color=palette,xlab=xlab,ylab=ylab,FUN=FUN)

}

##Missing values f�r die PC entfernen und nachher wieder zur�ck geben



zonalmean<-function(data,...)
{
	result<-list()
	result$zmean=latmean(data,...)
	result$lat = getlat(data)
                result$time=time(data)
	return(result)
}

