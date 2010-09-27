#climate indices (17.August.06, tlaepple@awi-bremerhaven.de)
#input time series field (pField)
#output time series (pTs) 

#NINO3 in SST (K) anomalies
index.nino3 <- function(sst)
{
  data<-selspace(sst,lon1=360-150,lon2=360-90,lat1=-5,lat2=5)
  result<-applyspace(scale(data,scale=FALSE),mean)
  attr(result,"name")<-paste("NINO 3 index",getname(sst))
  return(result)
}


#NINO3.4 in SST (K) anomalies
index.nino3.4 <- function(sst)
{
  data<-selspace(sst,lon1=360-170,lon2=360-120,lat1=-5,lat2=5)
  result<-applyspace(scale(data,scale=FALSE),mean)
  attr(result,"name")<-paste("NINO 3.4 index",getname(sst))
  return(result)

}




#NINO3.4 in SST (K) anomalies
index.nino1.2 <- function(sst)
{
  data<-selspace(sst,lon1=360-90,lon2=360-80,lat1=-10,lat2=0)
  result<-applyspace(scale(data,scale=FALSE),mean)
  attr(result,"name")<-paste("NINO 1+2 index",getname(sst))
  return(result)

}



#NINO3.4 in SST (K) anomalies
index.nino4 <- function(sst)
{
  data<-selspace(sst,lon1=360-160,lon2=360-150,lat1=-5,lat2=5)
  result<-applyspace(scale(data,scale=FALSE),mean)
  attr(result,"name")<-paste("NINO 4 index",getname(sst))
  return(result)

}

index.nino.tni<-function(sst)
{
  result<-scale(index.nino1.2(sst)-index.nino4(sst))
  attr(result,"name")<-paste("TNI index",getname(sst))
  return(result)

}


#NAO Index EOF-based
index.nao<-function(slp,plot=FALSE,pattern=FALSE)
{
	slpdata.naoregion<-selspace(slp,lon1=270,lon2=40,lat1=20,lat2=80) #select atlantic region
	result<-prcompO.pField(scale(detrend(slpdata.naoregion)))
	if (selspace(result$eof[1,],lat1=35.2,lon1=340.2) < 0)
		{
			warning("Sign of NAO was adapted to get a positive EOF over the Azores")
			result$pc[,1]<-result$pc[,1]*(-1)
		}
	if (plot) plot(result$eof[1,],"NAO / EOF1")
	attr(result$pc,"name")<-paste("NAO index/PC1",getname(slp))
	if (pattern) return(list(ts=scale(result$pc[,1]),pattern=result$eof[1,]))
			else return(scale(result$pc[,1]))

}

#original Tahiti-Darwin
index.soi <- function(slp)
{
  tahiti<-selspace(slp,lat1=(-17),lon1=(360-149))
  #tahiti<-selspace(slp,lat1=(-17),lon1=(360-149))

  darwin<-selspace(slp,lat1=(-12),lon1=(0+131))

  result<-scale(tahiti)-scale(darwin)
  attr(result,"name")<-paste("SOI index",getname(slp))
  return(result)

}

#adapted positions to get max correlation in the SONNE runs
index.soimodel <- function(slp)
{
  tahiti<-selspace(slp,lat1=0,lon1=239)
  darwin<-selspace(slp,lat1=-2,lon1=125)

  result<-scale(tahiti)-scale(darwin)
  attr(result,"name")<-paste("SOI index",getname(slp))
  return(result)

}

#PC1 definition

index.ao <- function(slp)
{
	slpdata.aoregion<-selspace(slp,lon1=0,lon2=360,lat1=20,lat2=90) #select atlantic region
	result<-prcompO.pField(scale(detrend(slpdata.aoregion)))
	if (selspace(result$eof[1,],lat1=35,lon1=340) < 0)
		{
			warning("Sign of AO was adapted to get a positive EOF over the Azores")
			result$pc[,1]<-result$pc[,1]*(-1)
		}
	
	result<-result$pc[,1]
	attr(result,"name")<-paste("AO index/PC1",getname(slp))
	return(result)

}

#Z(20N,160W) - Z(45N,165W) + Z(55N,115W) - Z(30N,85W) Wallace Definition

index.pna <- function(gph)
{
#Z(20N,160W) - Z(45N,165W) + Z(55N,115W) - Z(30N,85W)

   p1 <-selspace(gph,lat1=(20),lon1=(360-160))
   p2 <-selspace(gph,lat1=(45),lon1=(360-165))

   p3 <-selspace(gph,lat1=(55),lon1=(360-115))
   p4 <-selspace(gph,lat1=(30),lon1=(360-85))

 
  result<-0.25*(scale(p1)-scale(p2)+scale(p3)-scale(p4))
  attr(result,"name")<-paste("PNA index",getname(gph))
  return(result)

}




