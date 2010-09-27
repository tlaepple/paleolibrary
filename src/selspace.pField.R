#Selspace, 04.July.06  tlaepple@awi-bremerhaven.de


#boundaries are still missing
selspace <- function(data,lat1,lat2=NULL,lon1,lon2=NULL,tolLat=NULL,tolLon=NULL)
{
  
  temp<-attributes(data)
  stime<-time(data)

if (sum(c(is.null(lat2),is.null(lon2))) == 2)
  {
    #one point is choosen -> return time series on a specific point
    #derived by bilinear interpolation 



    #if no tolerance values are supplied take 1/2 box as tolerance
    if (is.null(tolLon)) tolLon<-abs((temp$lon[2]-temp$lon[1])/2)
    if (is.null(tolLat)) tolLat<-abs((temp$lat[2]-temp$lat[1])/2)


    #First search the surrounding longitudes and latitudes
    #order the latitudes and longitudes in a rising order....

    #attention... midpoints are given...

    indexLat2<-which(temp$lat>=lat1)[1]
    indexLon2<-which(temp$lon>=lon1)[1]
    indexLat1<-rev(which(temp$lat<=lat1))[1]
    indexLon1<-rev(which(temp$lon<=lon1))[1]

    #four special cases... 
    if ((is.na(indexLat1)+is.na(indexLat2)+is.na(indexLon1)+is.na(indexLon2))>0) stop("Point outside boundaries... this feature is not implemented yet")

    if (is.na(indexLat1))
	{
		if ((temp$lat[1]-lat1)<=tolLat)
		{warning("Point outside boundaries... but inside tolerance")
             indexLat1<-indexLat2
		}
			else stop("Point outside boundaries, increase tolerance")
	}
 
	if (is.na(indexLat2))
	{
		if ((lat1-temp$lat[length(temp$lat)])<=tolLat)
		{warning("Point outside boundaries... but inside tolerance")
             indexLat2<-indexLat1
		}
			else stop("Point outside boundaries, increase tolerance")
	}
#NOT FINISHED
    if (is.na(indexLon1)) #left boundary
	{
            #Test against 
		if ((temp$lat[1]-lat1)<tolerance)
		{warning("Point outside boundaries... but inside tolerance")
             indexLat1<-indexLat2
		}
			else stop("Point outside boundaries, increase tolerance")
	}
#NOT FINISHED
    if (is.na(indexLon2))
	{
		if ((lat1-temp$lat[0])<tolerance)
		{warning("Point outside boundaries... but inside tolerance")
             indexLat1<-indexLat2
		}
			else stop("Point outside boundaries, increase tolerance")
	}
nLon<-length(temp$lon)

#here the interpolation starts
d11<-data[,c2t1(indexLat1,indexLon1,nLon)]
d12<-data[,c2t1(indexLat2,indexLon1,nLon)]
d21<-data[,c2t1(indexLat1,indexLon2,nLon)]
d22<-data[,c2t1(indexLat2,indexLon2,nLon)]

ex=  (lon1-temp$lon[indexLon1])/(temp$lon[indexLon2]-temp$lon[indexLon1])
ey=  (lat1-temp$lat[indexLat1])/(temp$lat[indexLat2]-temp$lat[indexLat1])
#print(paste("DEBUG",indexLat1,indexLon1))


if ((!is.finite(ex))&(!is.finite(ey))) intpoldata<-d11 #we are on the point, no interpolation
else
if (!is.finite(ex))
{
	 intpoldata<-d11+(d12-d11)*ey #only latitudonal interpolation
}
   else 
{  
	if (!is.finite(ey)) 
	{
		intpoldata<-d11+(d21-d11)*ex #only longitudonal interpolation
	}
	else intpoldata <- (1-ex)*(1-ey)*d11 + (1- ex)*(ey)*d12 + ( ex)*(1-ey)*d21 + (ex*ey)*d22  
}

#create time series
result<-pTs(intpoldata,time(data),lat1,lon1,getname(data),gethistory(data),date=FALSE)

hist<-paste("selspace: lat=",lat1," lon=",lon1,sep="")

return(addhistory(result,hist));            

  }
  
else
  if  (sum(c(is.null(lat2),is.null(lon2))) == 0)
    {
#area is choosen
#choose all
  if (lat1 > lat2)
    indexLat<-(((temp$lat>=lat1) & (temp$lat <=90))
             |((temp$lat<=lat2) & (temp$lat >=-90)))
else
  indexLat<-((temp$lat>=lat1) & (temp$lat <=lat2))

if (lon1 > lon2)
    indexLon<-(((temp$lon>=lon1) & (temp$lon <=360))
             |((temp$lon<=lon2) & (temp$lon >=0)))
else
  indexLon<-((temp$lon>=lon1) & (temp$lon <=lon2))

  
  nLat<-length(temp$lat)
  nLon<-length(temp$lon)

  index2D<-NULL
  for (iLat in 1:nLat) 
      {
        if (indexLat[iLat]) index2D<-c(index2D,indexLon)
        else index2D<-c(index2D,rep(FALSE,nLon))
         }

  class(data)<-NULL
  subset<-data[,index2D]
  result<-pField(t(subset),stime,temp$lat[indexLat],temp$lon[indexLon],name=getname(data),gethistory(data),date=FALSE)
  hist<-paste("selspace: lat=",lat1,":",lat2," lon=",lon1,":",lon2,sep="")
  return(addhistory(result,hist));           

}
else stop("Either supply both edges or only one point ")
}
