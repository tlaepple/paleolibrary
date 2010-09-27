#concept for running calculations
#running correlation (+ significance), running variance

#input: ts1, 
roll.1 <- function(ts1,width,FUN,by=1,name=NULL,detrend=TRUE,scale=TRUE,... )
{
	nTimeGuess<-(end(ts1)[1]-start(ts1)[1]) / by
	StartTime<-(0:nTimeGuess)*by+start(ts1)[1]
	StartTime<-StartTime[StartTime<=end(ts1)[1]-width]
	nTime<-length(StartTime)
	
	newName<-paste(name,getname(ts1))
	result<-pTs(NULL,StartTime+(width/2),9999,9999,newName,gethistory(ts1),date=FALSE)
	
	for (i in 1:nTime) 
	{
		data<-window(ts1,StartTime[i],StartTime[i]+width)
		if (detrend) data<-detrend(data)
		if (scale) data<-scale(data)
		result[i]<-FUN(data, ...)
	}

	return(addhistory(result,newName))
}
roll.2 <- function(ts1,ts2,width,FUN,...,by=1,name=NULL,detrend=TRUE,scale=TRUE)
{
	#bring on same timebase
	start<-max(start(ts1)[1],start(ts2)[1])
	end<-min(end(ts1)[1],end(ts2)[1])
	print(paste("Common time period: ",start,end))
	ts1<-window(ts1,start,end)
	ts2<-window(ts2,start,end)

	nTimeGuess<-(end-start) / by
	StartTime<-(0:nTimeGuess)*by+start
	StartTime<-StartTime[StartTime<=end-width]
	nTime<-length(StartTime)
	
	newName<-paste(name,getname(ts1),getname(ts2))
	result<-pTs(NULL,StartTime+(width/2),9999,9999,newName,gethistory(ts1),date=FALSE)
	
	
	for (i in 1:nTime) 
	{
		data1<-window(ts1,StartTime[i],StartTime[i]+width)
		data2<-window(ts2,StartTime[i],StartTime[i]+width)
		if (detrend) {data1<-detrend(data1);data2<-detrend(data2)}
		if (scale) {data1<-scale(data1);data2<-scale(data2)}
		result[i]<-FUN(data1,data2, ...)
	}


	return(addhistory(result,newName))
}

cor.sig <- function(ts1,ts2,pval=0.05)
{
	temp<-cor.test(ts1,ts2)
	if (temp$p.value > pval) return(NA)
		else return(temp$estimate)
	
}





rollmean.k<-function(x,k,na.pad = FALSE, ...)
{
  newindex<-seq(floor(k/2)+1,by=k,length.out=length(x)/k)
  newtime<-time(x)[newindex]
  index<-!is.na(newtime)
  result<-pTs(NA,newtime[index],getlat(x),getlon(x),getname(x))
  result[]<-rollmean.default(x,k,na.pad=TRUE, ...)[newindex[index]]
  return(na.omit(result))

}


