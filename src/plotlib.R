 #plotting library, 15.August.06  tlaepple@awi-bremerhaven.de
#Update 15.04.08
#16.04.2008: axes=axes

#plot a eof calculated from multiple timeseries
#ts for getting the lat's and lon's
ploteof <- function(pcomp,ts,i=1,minLat=NULL,maxLat=NULL,minLon=NULL,maxLon=NULL)
{

 if (is.null(c(minLat,maxLat,minLon,maxLon)))
 {
	minLat<-min(getlat(ts))
	maxLat<-max(getlat(ts))
	extraLat<-(maxLat-minLat)/8
	minLat<-minLat-extraLat
	maxLat<-maxLat+extraLat

	minLon<-min(getlon(ts))
	maxLon<-max(getlon(ts))

	extraLon<-(maxLon-minLon)/8
	minLon<-minLon-extraLon
	maxLon<-maxLon+extraLon
 }

	#Give colors to the points
	a<-vector()
	a[pcomp$eof[,i]>0]<-"red"
	a[pcomp$eof[,i]<0]<-"blue"
	titlevar<-sprintf("%2.1f",pcomp$var[i]/pcomp$varsum*100)
	title=paste("EOF ",i,"  Expl.Var:",titlevar,"%",  sep="")
	plot(getlon(ts),getlat(ts),ylim=c(minLat,maxLat),xlim=c(minLon,maxLon),lwd=4,xlab="EAST",ylab="NORTH",col=a,main=title)
	text(getlon(ts),getlat(ts),sprintf("%1.2f",pcomp$eof[,i]),pos=1,offset=c(1,1))

	text(getlon(ts),getlat(ts),getname(ts),pos=3,offset=c(1,1))

	addland(col="black")
	grid()

}

plotwind<-function(data,title=NULL,zlim=range(data,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=rbow,FUN=NULL, ...)
  {
    tmp<-attributes(data)
    dim(data)<-c(length(getlat(data)),length(getlon(data)))
    sSub<-""

filled.contour(tmp$lat,tmp$lon,data,zlim=zlim,nlevels=nlevels,levels=levels,color=palette,plot.title={
        title(main=title,sub=sSub);
        if (!is.null(FUN)) FUN(tmp$lat,tmp$lon,data)
        addland(col="black");
        grid()})
      }



#Plots an index
plotindex <- function(data,rintervall=5,lwd=3,lwd1=2,main=NULL)
{
	if (is.null(main)) main<-getname(data)
	data.plus<-data
	data.minus<-data
	data.plus[data<=0]<-0
	data.minus[data>=0]<-0
	plot.ts(data,type="n",xlab="time",ylab="std",main=main)
	lines(data.minus,type="h",col="blue",lwd=lwd)
	lines(data.plus,type="h",col="red",lwd=lwd)
	lines(rollmean(data,5),col="black",lwd=lwd1)
}



plot.pTs<-function(x, plot.type = c("multiple", "single"), xy.labels,
    xy.lines, panel = lines, nc, yax.flip = FALSE, mar.multi = c(0,
        5.1, 0, if (yax.flip) 5.1 else 2.1), oma.multi = c(6,
        0, 5, 0), axes = TRUE, main=NULL,...)
{

 if (!is.null(ncol(x))) colnames(x)<-strtrim(getname(x),8) else if (is.null(main)) main=getname(x)
 plot.ts(x, NULL, plot.type, xy.labels,
    xy.lines, panel, nc, yax.flip , mar.multi, oma.multi,axes = axes,main=main, ...)
}

plot.pField <- function(x, ...) plotmap(x, ...)
plotmap<-function(plotdata, ...) UseMethod("plotmap")

myfun1<-function(sTitle,sSub)
{
       title(main=sTitle,sub=sSub)
       addland(col="black")
       grid()
}

myfun2<-function(sTitle,sSub,lat,lon,plotdata)
{ title(main=sTitle,sub=sSub)
contour(lon,lat,plotdata,col="white",zlim=c(min(plotdata),0),lty=1,lwd=2,nlevels=5,add=TRUE);
contour(lon,lat,plotdata,col="grey20",zlim=c(0,max(plotdata)),lty=1,lwd=2,nlevels=5,add=TRUE);
addland(col="black")
grid()
}

#checks the data
#wrap the data to get a continous lat/lon field
#reverse latitudes if needed
plot.preparation <- function(plotdata,shift=F,long=F)
{
        temp<-attributes(plotdata)
  	if (prod(dim(plotdata)) != length(temp$lon)*length(temp$lat)) stop("N(data) != N(lat)*N(lon)")

        plotdata<-matrix(plotdata,length(temp$lon),length(temp$lat)) #make a 2D array

        #arrange to get a continous field
        d<-diff(temp$lon)
        if (max(d) > (min(d)+0.01))
          { nlon<-length(temp$lon)
            edgelon<-which(d==max(d))
            plotdata<-rbind(plotdata[(edgelon+1):nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon+1):nlon],temp$lon[1:edgelon]+360)
          }

	if (shift) {
		nlon<-length(temp$lon)
            edgelon<-which(temp$lon>180)[1]
            plotdata<-rbind(plotdata[(edgelon-1):edgelon,],plotdata[(edgelon+1):nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon-1):edgelon]-360,temp$lon[(edgelon+1):nlon]-360,temp$lon[1:edgelon])
	}


	if (long)
	{
		nlon<-length(temp$lon)
            edgelon<-which(temp$lon>180)[1]
            plotdata<-rbind(plotdata[1:nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[1:nlon]-360,temp$lon[1:edgelon])


	}

        if (temp$lat[2] < temp$lat[1])    #if the latitudes are from + to -, reverse them
         {
          temp$lat<-rev(temp$lat)
          plotdata<-plotdata[,rev(seq(len=ncol(plotdata)))]
        }
        return(list(data=plotdata,lat=temp$lat,lon=temp$lon))
}

icecontour<-function(lon,lat,data)
  {
    contour(lon,lat,data,levels=0.5,add=TRUE,drawlabels=FALSE,col="blue",lwd=2)
  }
normcontour1<-function(lon,lat,data)
  {

    levels<-(1:5)/5
    contour(lon,lat,data,col="black",zlim=c(min(data),0),lty=2,lwd=1,levels=-1*levels,add=TRUE,drawlabels=FALSE)
    contour(lon,lat,data,col="black",zlim=c(0,max(data)),lty=1,lwd=1,levels=levels,add=TRUE,drawlabels=FALSE)
    #contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=FALSE,col="black",lwd=2)
  }

slpcontour<-function(lon,lat,data)
  {


     contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=T,col="black",lwd=1)
  }



slp.diff.contour<-function(lon,lat,data)
  {
     levelsp=c(2,4,6,8,10,12,14,16)
    levelsm= -1*levelsp

    contour(lon,lat,data,col="black",lty=1,lwd=1,levels=levelsp,add=TRUE,drawlabels=T)
    contour(lon,lat,data,col="black",lty=2,lwd=1,levels=levelsm,add=TRUE,drawlabels=T)
    contour(lon,lat,data,col="grey",lty=1,lwd=2,levels=c(0),add=TRUE,drawlabels=T)



  }



normcontour0<-function(lon,lat,data)
  {

    levels<-(1:8)/80
    contour(lon,lat,data,col="black",zlim=c(min(data),0),lty=2,lwd=1,levels=-1*levels,add=TRUE,drawlabels=FALSE)
    contour(lon,lat,data,col="black",zlim=c(0,max(data)),lty=1,lwd=1,levels=levels,add=TRUE,drawlabels=FALSE)
    #contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=FALSE,col="black",lwd=2)
  }

pnacontour<-function(lon,lat,data)
  {

    levelsp=c(0.25,0.5,0.75)
    levelsm=c(-0.25,-0.5,-0.75)

    contour(lon,lat,data,col="black",lty=2,lwd=1,levels=levelsp,add=TRUE,drawlabels=FALSE)
    contour(lon,lat,data,col="black",lty=1,lwd=1,levels=levelsm,add=TRUE,drawlabels=FALSE)
    contour(lon,lat,data,col="grey",lty=1,lwd=2,levels=c(0),add=TRUE,drawlabels=FALSE)

    #contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=FALSE,col="black",lwd=2)
  }


corcontour<-function(lon,lat,data)
  {

    levelsp=c(0.2,0.4,0.6)
    levelsm=c(-0.2,-0.4,-0.6)

    contour(lon,lat,data,col="black",lty=2,lwd=1,levels=levelsp,add=TRUE,drawlabels=TRUE)
    contour(lon,lat,data,col="black",lty=1,lwd=1,levels=levelsm,add=TRUE,drawlabels=TRUE)
    contour(lon,lat,data,col="grey",lty=1,lwd=2,levels=c(0),add=TRUE,drawlabels=TRUE)

    #contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=FALSE,col="black",lwd=2)
  }



gphcontour<-function(lon,lat,data)
  {

    levelsp=c(5,10,15,20)
    levelsm=c(-5,-10,-15,-20)

    contour(lon,lat,data,col="black",lty=2,lwd=1,levels=levelsp,add=TRUE,drawlabels=TRUE)
    contour(lon,lat,data,col="black",lty=1,lwd=1,levels=levelsm,add=TRUE,drawlabels=TRUE)
    contour(lon,lat,data,col="grey",lty=1,lwd=2,levels=c(0),add=TRUE,drawlabels=TRUE)

    #contour(lon,lat,data,nlevels=10,add=TRUE,drawlabels=FALSE,col="black",lwd=2)
  }


cwind<-function(lon,lat,data)
  {

    levelsp<-(1:5)*8
    levelsm<--2*(1:10)

    contour(lon,lat,data,col="black",zlim=c(min(data),0),lty=2,lwd=1,levels=levelsm,add=TRUE,drawlabels=FALSE)
    contour(lon,lat,data,col="black",zlim=c(0,max(data)),lty=1,lwd=1,levels=levelsp,add=TRUE,drawlabels=FALSE)

  }

rbow <- function (n, s = 1, v = 1, start = 0, end = 0.7, gamma = 1)
{
    if ((n <- as.integer(n[1])) > 0) {
        if (start == end || any(c(start, end) < 0) || any(c(start,
            end) > 1))
            stop("'start' and 'end' must be distinct and in [0, 1].")
        result<-hsv(h = seq(start, ifelse(start > end, 1, 0) + end, length = n)%%1,
            s, v, gamma)
        reverse(result)
    }
    else character(0)
}


##########

plotmap.pField <- function(plotdata,main=NULL,zlim=range(plotdata,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=NULL,FUN=NULL,shift=F,long=F,xlim=NULL,stype=0,sSub="", ...)
{

	temp<-attributes(plotdata)

	if (is.null(sSub)) if (time(plotdata) != 9999) sSub<-paste("time:",format(time(plotdata)))
	if (is.null(main)) main<-temp$name
	gridcolor="lightgray"

	if (stype == 1) {
		shift=T
		xlim=c(-180,180)
		if (is.null(palette)) {
				palette=colorRampPalette(c("violetred4","blue","steelblue", "lightgreen","white", "yellow","orange","red","brown"))
				gridcolor="black"
				}

			}

	if (stype == 2) {

		if (is.null(palette)) {
				palette=colorRampPalette(c("violetred4","blue","steelblue", "lightgreen","white", "yellow","orange","red","brown"))
				gridcolor="black"
				}

			}


	if (is.null(palette)) { palette=rbow;}

	#Changes from the 1D vector to a 2D field
      tmp<-plot.preparation(plotdata,shift,long)

	#Avoid longitudes > 360
	if (max(tmp$lon)>360) tmp$lon<-tmp$lon-360

	if (is.null(xlim))  xlim = range(tmp$lon, finite = TRUE)

	if (stype == 1)
	{
filled.contour.own(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,xlim=xlim,color=palette,plot.title={
        title(main=main,sub=sSub);
	addland(col="black");
        if (!is.null(FUN)) FUN(tmp$lon,tmp$lat,tmp$data)

        },plot.axes=axes.stype(gridcolor,tmp$lat,tmp$lon))
      }

	else
	filled.contour.own(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,xlim=xlim,color=palette,plot.title={
        title(main=main,sub=sSub);
	addland(col="black"); grid()
        if (!is.null(FUN)) FUN(tmp$lon,tmp$lat,tmp$data)

       })
      }


addland<-function (col = "grey50", lwd = 1)
{
    data("addland1")
    lines(lon.cont, lat.cont, type = "l", col = col, lwd = lwd)
    lon.cont <- lon.cont + 360
    lines(lon.cont, lat.cont, type = "l", col = col, lwd = lwd)
	lon.cont <- lon.cont - 720
    lines(lon.cont, lat.cont, type = "l", col = col, lwd = lwd)

}

axes.stype<-function(gridcolor,lat,lon)
{

if (min(lat)<0) {
	labels.lat<-c("60S","30S","EQ","30N","60N")
	at.lat<-c(-60,-30,0,30,60)
	}
else
{
	labels.lat<-c("10S","30N","50N","70N")
	at.lat<-c(10,30,50,70)
}
labels.lon<-c("180","90W","0","90E","180")
at.lon=c(-180,-90,0,90,180)
title(main = "", xlab = "", ylab = "")
            Axis(lat, side = 1,at=at.lon,labels=labels.lon)
            Axis(lon, side = 2,at=at.lat,labels=labels.lat)
abline(h=at.lat,col = gridcolor, lty = "dotted")
abline(v=at.lon,col = gridcolor, lty = "dotted")

}

plotmap.pFieldb <- function(plotdata,plotdata2,title=NULL,zlim=range(plotdata,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=rbow,FUN=NULL, ...)
{
	temp<-attributes(plotdata)
      sSub<-NULL
	if (time(plotdata) != 9999) sSub<-paste("time:",format(time(plotdata)))
	if (is.null(title)) title<-temp$name
        tmp<-plot.preparation(plotdata)
	  tmp2<-plot.preparation(plotdata2)

	filled.contour(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,color=palette,plot.title={
        title(main=title,sub=sSub);
        if (!is.null(FUN)) FUN(tmp2$lon,tmp2$lat,tmp2$data)
        addland(col="black");
        grid()})
      }


plotmapc.pField <- function(plotdata,sTitle=NULL, ...)
{
	temp<-attributes(plotdata)
 	sSub<-paste("time:",format(time(plotdata)))
	if (is.null(sTitle)) sTitle<-temp$name

	if (prod(dim(plotdata)) != length(temp$lon)*length(temp$lat)) stop("N(data) != N(lat)*N(lon)")

        plotdata<-matrix(plotdata,length(temp$lon),length(temp$lat)) #make a 2D array

        #arrange to get a continous field
        d<-diff(temp$lon)
        if (max(d) > (min(d)+0.01))
          { nlon<-length(temp$lon)
            edgelon<-which(d==max(d))
            plotdata<-rbind(plotdata[(edgelon+1):nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon+1):nlon],temp$lon[1:edgelon]+360)
          }


        if (temp$lat[2] < temp$lat[1])    #if the latitudes are from + to -, reverse them
         {
          temp$lat<-rev(temp$lat)
          plotdata<-plotdata[,rev(seq(len=ncol(plotdata)))]
                    }

	 levels=(0:6)*0.01
	 contour(temp$lon,temp$lat,plotdata,col="blue",zlim=c(min(plotdata),0),lty=1,lwd=2,levels=-1*levels)
         contour(temp$lon,temp$lat,plotdata,col="red",zlim=c(0,max(plotdata)),lty=2,lwd=2,levels=levels,add=TRUE)
     title(main=sTitle,sub=sSub);addland(col="black");grid();
      }

plotcont.pField <- function(plotdata,sTitle=NULL, ...)
{
	temp<-attributes(plotdata)
 	sSub<-paste("time:",format(time(plotdata)))
	if (is.null(sTitle)) sTitle<-temp$name

	if (prod(dim(plotdata)) != length(temp$lon)*length(temp$lat)) stop("N(data) != N(lat)*N(lon)")

        plotdata<-matrix(plotdata,length(temp$lon),length(temp$lat)) #make a 2D array

        #arrange to get a continous field
        d<-diff(temp$lon)
        if (max(d) > (min(d)+0.01))
          { nlon<-length(temp$lon)
            edgelon<-which(d==max(d))
            plotdata<-rbind(plotdata[(edgelon+1):nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon+1):nlon],temp$lon[1:edgelon]+360)
          }


        if (temp$lat[2] < temp$lat[1])    #if the latitudes are from + to -, reverse them
         {
          temp$lat<-rev(temp$lat)
          plotdata<-plotdata[,rev(seq(len=ncol(plotdata)))]
                    }

      	 contour(temp$lon,temp$lat,plotdata,col="blue",zlim=c(min(plotdata),0),lty=1,lwd=2,nlevels=5)
         contour(temp$lon,temp$lat,plotdata,col="red",zlim=c(0,max(plotdata)),lty=2,lwd=2,nlevels=5,add=TRUE)
     title(main=sTitle,sub=sSub);addland(col="black");grid();
      }


makefilm<-function(data,startdate,enddate,...,avrg=11,step=5,prefix="ani_",anomaly=FALSE)
{

meandata<-applyspace(data,mean)
meandata<-na.omit(rollmean(meandata,avrg))
ypos<-max(meandata)

if (anomaly) data<-data-applytime(data,mean)

nStep<-(enddate-startdate) %/% step
starttime=(1:nStep)*step + startdate
starttime<-starttime[starttime<(enddate-avrg)]


i<-100
for (it in starttime)
{
  i<-i+1
  tdata<-window(data,it,it+avrg)
  result<-applytime(tdata,mean)
  jpeg(paste(prefix,i,".jpg",sep=""),width=600,height=600,quality=90)

  plot(result,title=paste(it,"-",it+avrg),...)

  dev.off()

}
}


filled.contour.own<-function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)),
    z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE),
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) -
        1), plot.title, plot.axes, key.title, key.axes, asp = NA,
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes,
    ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
        yaxs = "i")
  #  rect(0, levels[-length(levels)], 1, levels[-1], col = col)

    delta<-(levels[length(levels)]-levels[1])/(length(levels)-1)
    breaks<-delta*(0:(length(levels)-1))+levels[1]
    rect(0, breaks[-length(levels)], 1, breaks[-1], col = col)

    if (missing(key.axes)) {
        if (axes)
		{    #use modified axes
            	axis(4,labels=levels,at=breaks)
		}
    }
    else key.axes
    box()
    if (!missing(key.title))
        key.title
    mar <- mar.orig
    mar[4] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
        stop("no proper 'z' matrix specified")
    if (!is.double(z))
        storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels),
        col = col))
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}
 rbow.symetric=colorRampPalette(c("violetred4","blue","steelblue", "lightgreen","white", "yellow","orange","red","brown"))

