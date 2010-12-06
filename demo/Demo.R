 # Aim: Demonstration of the PaleoLibrary
# Date:
# Documentation:

#Changes:

path="e:/data/paleoLibrary/src/"
source(paste(path,"header.R",sep=""))


tair<-read_ncep("e:/data/pool/ncep/air.annual.nc",varname="air")

slp<-read_ncep("e:/data/pool/ncep/slp.ann.nc",varname="slp")

result<-read_data("e:/data/pool/gridded/annual/air.ann.nc",varname="air")

plot(tair[,1])
time(tair)

index.nao(slp,plot=TRUE)
plotindex(index.nao(slp))

mean(window(applyspace(tair,mean),1948,1960))

mean(window(applyspace(tair,mean),2003,2007))

ipcc.mean<-read_ncep("e:/data/HolFut/ens.ann.nc",varname="tas")
ipcc.mean[ipcc.mean>1e18]<-273.15

clim.ipcc<-applytime(window(ipcc.mean,1930,1960),mean)
clim.ipcc.x<-pField(rep(clim.ipcc,dim(ipcc.mean)[1]),time(ipcc.mean),getlat(ipcc.mean),getlon(ipcc.mean))

karte(ipcc.mean,2080,2090)
ipcc.mean.anom<-ipcc.mean-clim.ipcc.x

mean(clim.ipcc)
mean(

     karte(ipcc.mean.anom,2020,2040,dlim=c(-10,10))
     karte(lgm,250,280,dlim=c(-20,20))

r<-read_ncep("e:/data/HolFut/air.noaa.ann.nc",varname="air")
lgm<-read_ncep("e:/data/ARCH/coupled/LGM_985/temp.ann.nc",varname="temp2")
pi<-read_ncep("e:/data/ARCH/coupled/PI_986/temp.ann.nc",varname="temp2")

clim.pi<-applytime(pi,mean)
clim.pi.300<-pField(rep(clim.pi,300),time(lgm),getlat(lgm),getlon(lgm))

lgm<-window(lgm,250,300)

lgm<-(lgm-clim.pi.300)

karte(lgm,200,250,dlim=c(-20,20))

zeitserie(lgm,-90,0,90,360)
zeitserie(clim.pi.300,-90,0,90,360)




plot(r[190,]-r[1,])
r[r>1000]<-NA
plot(applytime(r,mean))

plot(tair[,1])
par(mfcol=c(2,2))

plotindex(applyspace(tair,mean))
abline(v=c(1992,1981))

clim<-applytime(window(tair,1950,1970),mean)
plot(window(tair,1998,1998)-clim,zlim=c(-1,2))

tbremen<-selspace(r,lat1=53,lon1=10)
plot(tbremen)

plot(tbremen,ylab="Temperature, Grad Celsisus",xlab="Jahr")
cbremen<-cor.pTs(detrend(tbremen),tair)
plot(cbremen,stype=1,zlim=c(-1,1))

r<-prcompO.pField(tair)


karte<-function(data,startjahr,endjahr,dlim=range(data))
{
    temp<-applytime(window(data,startjahr,endjahr),mean)
    temp<-selspace(temp,lat1=-86,lat2=86,lon1=0,lon2=360)
    plot(temp,stype=1,zlim=dlim,main=paste(startjahr,"-",endjahr))
}


zeitserie<-function(data,lat1,lon1,lat2=NULL,lon2=NULL,drange=range(data),add=FALSE)
{
    lon1t<-lon1
    lon2t<-lon2
    if (lon1 < 0) lon1t=lon1+360
  if (!is.null(lon2))  if (lon2 < 0) lon2t=lon2+360
    temp<-selspace(data,lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2)
    if (!is.null(dim(temp))) temp<-applyspace(temp,mean)

      if (!is.null(lon2)) main=paste(lat1,"-",lat2,"N   ",lon1,"-",lon2,"W",sep="") else main=paste(lat1,"N   ",lon1,"W",sep="")
    if (add) lines(temp,col="red",lwd=2)

       else  plot(temp,ylab="Temperaturanomalie (K)",xlab="Zeit (Jahr)",main=main,lwd=2)
}



zeitserie(tair,50,20,60,40)
zeitserie(tair,60,20,add=T)
zeitserie(tair,-90,0,90,360)
lines(applyspace(tair,mean),col="red")
karte(tair,1960,1970)
