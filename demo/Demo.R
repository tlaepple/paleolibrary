# Aim: Demonstration of the PaleoLibrary
# Date:
# Documentation: 

#Changes: 

path="e:/data/paleoLibrary/src/"
source(paste(path,"header.R",sep=""))




tair<-read_data("e:/data/pool/ncep/air.annual.nc")



plot(tair[,1])
tbremen<-selspace(tair,lat1=50,lon1=5)

plot(tbremen)
cbremen<-cor.pTs(detrend(tbremen),tair)
plot(cbremen,stype=1,zlim=c(-1,1))

r<-prcompO.pField(tair)

