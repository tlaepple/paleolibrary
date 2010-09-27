#proxytype library, 03.July.06 tlaepple@awi-bremerhaven.de
#17.august, correction in pTs type

#data[lat,lon,time]
pField <- function(data,time,lat=0,lon=0,name=" ",history=" ",date=TRUE)
{
#constants
TOL=1/400 #tolerance less than 1 day

#check data
      if (length(data) <=1 ) 
		if (is.null(data[1])) data<-rep(NA,length(lat)*length(lon)*length(time))
			else data<-rep(data[1],length(lat)*length(lon)*length(time))
			
	if (length(data) != (length(lat)*length(lon)*length(time))) 
		{	
			stop("nLat*nLon*nTime != N_Elements(data), if you want to create an empty field, supply data=NA")	
		}
	
#shape data in 2D array
	dim(data)<-c(length(lat)*length(lon),length(time))

      if (length(time) > 1) #real time series
      {
		if (abs(max(diff(time))-min(diff(time)))>TOL) stop("time steps are not equidistant")
      	result<-ts(t(data),start=time[1],deltat=(time[2]-time[1]))
	}
	else result<-ts(t(data),start=time[1]) #or only one time step

#put attributes and classes
	attr(result,'lat')<-lat
	attr(result,'lon')<-lon
	attr(result,'name')<-name
	if (date) attr(result,'history')<-paste(date(),history)
		else attr(result,'history')<-paste(history)


	attr(result,'oclass')<-class(result)
      attr(result,'nclass')<-c('pField','ts')
      class(result)<-attr(result,'nclass')
	
	invisible(result)
}





#data[x(t),i]

pTs<-function(data,time,lat=0,lon=0,name="",history="",date=TRUE)
{
#constants
TOL=1/400 #tolerance less than 1 day

 if (length(data) <=1 ) 
                if (is.null(data[1])) data<-matrix(NA,length(time),length(name))
                        else if (length(name) == 1) data<-rep(data[1],length(time)) else
					data<-matrix(data[1],length(time),length(name))

                
                
if (!is.null(ncol(data)) && (ncol(data)>1)) #multiple datasets
{       #check data
         if (nrow(data) != length(time)) stop("nTime != N_Elements(data)")

}
#check data
else if (length(data) != (length(time))) stop("nTime != N_Elements(data)")
        
#shape data in 2D array
        
      if (length(time) > 1) #real time series
      {
                if (abs(max(diff(time))-min(diff(time)))>TOL) stop("time steps are not equidistant")
        result<-ts(data,start=time[1],deltat=(time[2]-time[1]))
        }
        else result<-ts(data,start=time[1]) #or only one time step

#put attributes and classes
        attr(result,'lat')<-lat
        attr(result,'lon')<-lon
        attr(result,'name')<-name

                if (date) attr(result,'history')<-paste(date(),history)
                else attr(result,'history')<-paste(history)


        attr(result,'oclass')<-class(result)
      attr(result,'nclass')<-c('pTs','ts')
      class(result)<-attr(result,'nclass')
        
        invisible(result)
}
Ops.pField <- function (e1, e2) 
{
LIMIT<-1
TLIMIT<-1
#Field and scalar, scalar and Field... simply take the normal operator... and change history
#Field and Field: check if lat/lons are compatible... than transform in time series and back in pField
#.Generic gives the the name of the operator

nField <- 0
if (is_pField(e1)) { nField<-nField+1;name1<-attr(e1,"name")}
	else name1<-e1[1]

if (is_pField(e2)) {nField<-nField+1;name2<-attr(e2,"name")}
	else name2<-e2[1]

newname<-paste(name1,.Generic,name2)

if (nField == 2) 
{
	temp1<-attributes(e1)
	temp2<-attributes(e2)
	time1<-time(e1)

	if ((length(temp1$lat) != length(temp1$lat)) |(length(temp1$lon) != length(temp1$lon))) stop("grids not compatible")
	if ((sum(abs(temp1$lat-temp2$lat))+sum(abs(temp1$lon-temp2$lon))) > LIMIT) stop("grids not compatible")

	if (abs(sum(c(time(e1))-c(time(e2))))>TLIMIT) 	warning("Operator applied on two timebases, new timebase = first timebase")
				
#Trick, bring it on the same timebase... start as to be > 0 due to R-bug


	e1<-ts(e1,start=100)
	e2<-ts(e2,start=100)
	return(pField(t(NextMethod()),time1,temp1$lat,temp1$lon,newname,c(paste("A:",temp1$history),paste("B:",temp2$history),newname)))
}
else
{
e<-NextMethod(.Generic)
if (!is.null(attr(e,'nclass'))) class(e)<-attr(e,'nclass')
return(addhistory(e,newname))
}
}




Ops.pTs <- function (e1, e2) 
{
TLIMIT<-1
#Field and scalar, scalar and Field... simply take the normal operator... and change history
#Field and Field: check if lat/lons are compatible... than transform in time series and back in pField
#.Generic gives the the name of the operator

nField <- 0
if (is_pTs(e1)) { nField<-nField+1;name1<-attr(e1,"name")}
	else name1<-e1[1]

if (is_pTs(e2)) {nField<-nField+1;name2<-attr(e2,"name")}
	else name2<-e2[1]

newname<-paste(name1,.Generic,name2)

if (nField == 2) 
{
	temp1<-attributes(e1)
	temp2<-attributes(e2)
	time1<-time(e1)
#Trick, bring it on the same timebase... start as to be > 0 due to R-bug
	e1<-ts(e1,start=100)
	e2<-ts(e2,start=100)

	if (abs(sum(c(time(e1))-c(time(e2))))>TLIMIT) 	warning("Operator applied on two timebases, new timebase = first timebase")
	
	return(pTs(NextMethod(),time1,temp1$lat,temp1$lon,newname,c(paste("A:",temp1$history),paste("B:",temp2$history),newname)))
}
else
{
e<-NextMethod(.Generic)
if (!is.null(attr(e,'nclass'))) class(e)<-attr(e,'nclass')
return(addhistory(e,newname))
}
}







 
#if p1 and p2 [x,y] or only p [x] are given: give back value
      #if only p1 is given   [x,] :  pField: field at one time
      #if only p2 is given   [,x]: pTs   : univariate time series
   
	#über die Dimensionen arbeiten...
      #wenn [1:n,] dann length (dim=NULL) oder ncol (wenn dim != null) gleich ncol(x)
	#damit kann dies von [1:n] unterschieden werden 

      #Fälle zum testen: [a,b] [a:a1,b]  [a,b:b1] [a:a1,b:b1] 
      # [a,] [a:a1,] [,b] [,b:b1] [a] [a:a1] [date] .... das ganze auf normalem Feld.. und Feld mit einem


"[.pField"<-function(x,p1,p2,...)
{
  result<-NextMethod("[")	
  temp<-attributes(x)
  
  if (!missing(p2)) l2<-length(p2)

  if (!missing(p1))
  {
      l1<-length(p1)
	if (missing(p2)) 
	{     
		#check in [a] oder [a,]
		if ((is.null(dim(result)) & (length(result)==ncol(x))) | !is.null(dim(result))) 
		{  
		      #[a,]  give back field at times a:a1
			hist<-paste("[",p1[1],":",p1[l1],", ]",sep="")
		      if (length(p1) > 1) result<-t(result)
		      result<-pField(result,time(x)[p1],temp$lat,temp$lon,temp$name,temp$history,date=FALSE)
		}
		else
		{
			#[a]  directly give back value
			hist<-paste("[",p1[1],":",p1[l1],"]",sep="")						
		}
		
		
	}
  	else 
	{
	#[a,b] directly give back value or area
	hist<-paste("[",p1[1],":",p1[l1],",",p2[1],":",p2[l2],"]",sep="")
	}
 }
 else if (!missing(p2)) 
 {
	#[,b]  give back time series at position b
	hist<-paste("[,",p2[1],":",p2[l2],"]",sep="")

      pos<-c1t2(p2,length(temp$lon))
      result<-pTs(result,time(x),temp$lat[pos$lat],temp$lon[pos$lon],temp$name,temp$history,date=FALSE)
			
 }

  if (!is.null(attr(result,"history"))) return(addhistory(result,hist))  else return(result)
	

}


#[a,] time area a out of timeseries
#[a] time area a out of timeseries
#[a,b] point
#[,b] time series with index b

"[.pTs"<-function(x,p1,p2,...)
{
  result<-NextMethod("[")	
  temp<-attributes(x)
  
  if (!missing(p2)) l2<-length(p2)

  if (!missing(p1))
  {
      l1<-length(p1)
	if (missing(p2)) 
	{     
		#check in [a] oder [a,]
	
		      #[a,]  give back points at times a:a1
		       	hist<-paste("[",p1[1],":",p1[l1],",]",sep="")
                        result<-pTs(result,time(x)[p1],temp$lat,temp$lon,temp$name,temp$history,date=FALSE)

		     
			
		
	}
  	else 
	{
	#[a,b] directly give back value or area
	hist<-paste("[",p1[1],":",p1[l1],",",p2[1],":",p2[l2],"]",sep="")
	}
 }
 else if (!missing(p2)) 
 {
	#[,b]  give back time series at position b
	hist<-paste("[,",p2[1],":",p2[l2],"]",sep="")

      pos<-c1t2(p2,length(temp$lat))
      result<-pTs(result,time(x),temp$lat[p2],temp$lon[p2],temp$name[p2],date=FALSE)
			
 }
 
  if (!is.null(attr(result,"history"))) return(addhistory(result,hist))  else return(result)

}
