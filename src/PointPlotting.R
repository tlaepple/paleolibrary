# Aim: Demonstration code for plotting data on the maps
# Date: 31.July 09, Thomas Laepple
# Documentation: 

#Changes: 

path="e:/data/paleoLibrary/src/"
source(paste(path,"header.R",sep=""))

#Here I assume that temp.slope.echog is some test field already loaded (for example the field of the Holocene temperature slope of ECHO-G)


#Function for plotting points on maps using the same colors as inside the map
#lat,lon are vectors containing the latitude and longitude of the points
#value are the values of the points (e.g. slopes of the Foram data)

#pch = 19: filled circle,  pch=22 filled square
#lwd = size of the points
#Either provide the color levels (levels = ...)
#or provide the map data (data = ...) and it will calculate the levels form the map data

#By default, palette is the same one as used in the plot command with stype=2... 

addpoints<-function(lat,lon,value,pch=19,lwd=7,data=NULL,zlim=range(data,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=colorRampPalette(c("violetred4","blue","steelblue", "lightgreen","white", "yellow","orange","red","brown")))
{
	#if (is.null(levels)) stop("Levels of the plot are required")
	
	# First expand the levels to a higher resolution (here 10 sublevels for each)
	highres.levels<-(seq(levels)-1)*10+1
	allcolors<-palette(highres.levels[length(highres.levels)]) #Get all the colors
	colors<-allcolors[approx(levels,highres.levels,value)$y] #Get the colors of the points by linear interpolation
	points(lon,lat,col=colors,pch=pch,lwd=lwd)
}


#### Test code
                       

#Definition of some test datasets for the points and squares... 
latp<-c(-30,-20,-10,0,10,20)
lonp<-c(200,200,200,200,200,200)
valuep<-c(-2.5,-2,-1,0,-1.5,2)


latp1<-c(-30,-20,-10,0,10,20)
lonp1<-c(100,100,100,100,100,100)
valuep1<-c(-2.5,-2,-1,0,-1.5,2)


#Here define the colorlevels used in the map and the point plotting... you can adapt is as you want
levels=c(-3,-2,-1,-0.7,-0.5,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.5,0.7,1,2,3)


#Example for adding circles with the latitude latp, longitude lonp, and colors corresponding to the values valuep
#This is the normal plooting command, stype=2 for the violett-red colorscale, color levels are given by the levels vector defined above
#FUN = addpoints(...) is added to plot the points
#in addpoints, one has to give the latitudes,longitudes and values of the points to be added, the same color levels as in the main plot command (levels=levels)
#pch=19 means filled circles
plot(temp.slope.echog,stype=2,levels=levels,FUN=addpoints(latp,lonp,valuep,levels=levels,pch=19))

#Example for adding squares with the latitude latp1, longitude lonp1, and colors corresponding to the values valuep1
plot(temp.slope.echog,stype=2,levels=levels,FUN=addpoints(latp1,lonp1,valuep1,levels=levels,pch=22))

#For plotting multiple point arrays one can give multiple function to FUN by FUN = {function1;function2} 
plot(temp.slope.echog,stype=2,levels=levels,FUN={FUN=addpoints(latp,lonp,valuep,levels=levels,pch=19);addpoints(latp1,lonp1,valuep1,levels=levels,pch=22)})




#Example how it works without defining levels. Therefore you have to give the data to the addpoints function (here data = temp.slope.echog), that it knows the
#color range of the map
plot(temp.slope.echog,stype=2,FUN=addpoints(latp,lonp,valuep,data=temp.slope.echog,pch=19))
