#tlaepple, 17.August 06, ncol  bug

season.pTs <- function(ts,timewindow=c(1,11)/12,TOL=1/350,debug=FALSE)
{
	
	startyear<-floor(start.own(ts))
	endyear<-ceiling(end.own(ts))

	
	if (!is.null(ncol(ts)) && (ncol(ts)>1))
	{
	newTs<-pField(NULL,startyear:endyear,getlat(ts),getlon(ts),getname(ts),"season")


	for (i in startyear:endyear) 
	{
	

		if ((start.own(ts) <= (i+timewindow[1])) &&  (end.own(ts) >= (i+timewindow[2])) )
            {
		       newTs[i-startyear+1,]<-apply(window(ts,i+timewindow[1],i+timewindow[2]),2,mean)
		if (debug) {
				print(time(window(ts,i+timewindow[1],i+timewindow[2])))
				print(time(newTs)[i-startyear+1])
				}
		}

	}
	}
	else
	{newTs<-pTs(NULL,startyear:endyear,getlat(ts),getlon(ts),getname(ts),"season")

	for (i in startyear:endyear) 
	{
		if ((start.own(ts) <= (i+timewindow[1])) &&  (end.own(ts) >= (i+timewindow[2])) )
		{
                   newTs[i-startyear+1]<-mean(window(ts,i+timewindow[1],i+timewindow[2]))
			if (debug) {
				print(time(window(ts,i+timewindow[1],i+timewindow[2])))
				print(time(newTs)[i-startyear+1])
					}
		}
	}
	}
	#Auch die Nullen entfernen bei Feldern !
	return(na.omit(newTs))
#für alle mean rechnen
}

start.own <- function(ts) time(ts)[1]
end.own <- function(ts) time(ts)[length(time(ts))]

seascycle.pTs<-function(ts)
{
	start<-time(ts)[1]
	seascycle<-window(ts,start,start+0.99)
	nyears<-floor(time(ts)[length(ts)]-start)
      count<-1
      for (i in 1:nyears-1)
	{
		data<-window(ts,start+i,start+i+0.99)
		if (length(time(seascycle))== length(time(data))) 
		{
			seascycle<-seascycle+data
			count<-count+1
		}
	}
	return(seascycle/count)
}

#remove seasonal mean
rm_season.pTs<-function(ts)
{
	start<-time(ts)[1]
	seascycle<-window(ts,start,start+0.99)
	nyears<-floor(time(ts)[length(ts)]-start)
      count<-1
      for (i in 1:nyears-1)
	{
		data<-window(ts,start+i,start+i+0.99)
		if (length(time(seascycle))== length(time(data))) 
		{
			seascycle<-seascycle+data
			count<-count+1
		}
	}
	seascycle<-seascycle/count
	
	for (i in 1:nyears-1)
	{
		data<-window(ts,start+i,start+i+0.99)
		if (length(time(seascycle))== length(time(data))) 
		{
			window(ts,start+i,start+i+0.99)<-data-seascycle

		}
	}
return(ts)

}
