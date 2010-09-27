
timeser<-detrend(selspace(window(t2m,1500,1600),lat1=30,lon1=20))
result<-cortest.pTs(timeser,window(t2m,1500,1600))

plot(result[2,])
temp<-plot.preparation(result[2,],shift=T)
filled.contour(temp$data)
siglines<-contourLines(temp$lon,temp$lat,1-temp$data,levels=0.95) #Contour lines for Monsoon



plotSig<-function(sigline)
{
	col="grey60"
	for (i in 1:length(sigline)) polygon(sigline[[i]]$x,sigline[[i]]$y,angle=-30,density=30,col=col, border = col)
}



plot(result[1,],FUN=plotSig(siglines),stype=1,zlim=c(-1,1))



############################################

corcontur()

sigfield<-result[1,]
sigfield[result[2,]>0.05]<-NA

#Übergabe für shiftparameter in plotmap.pFieldb einbauen
#Contouren können in corcontour variiert werden

plotmap.pFieldb(sigfield,result[1,],FUN=corcontour,shift=T)  

#Beispiel:

levels=c(-30,-20,-10,-5,-3,-2,-1.5,-1,-0.5,-0.25,0.25,0.5,1,1.5,2,3,5,10,20,30)
plot((slopes[3,]-slopes[1,])/cfield.lin[2,],stype=1,levels=levels,main="Summer-Winterslope",FUN={addpoints();plotLagAndR2()})

