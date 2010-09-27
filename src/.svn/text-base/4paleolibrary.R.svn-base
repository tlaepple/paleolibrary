read_data<-function(FILENAME="C:/data/NCEP/DJFslp.nc",varname=NULL,name="",lonname=NULL,latname=NULL,missVal=c(-1e+20,1e+20)){
	temp.nc = open.ncdf(FILENAME)

	# read varname from temp.nc
	if(is.null(varname)){
		if(temp.nc$nvars>1){
			varnames<-c()
			for (i in 1:length(temp.nc$var)){
					varnames<-c(varnames,temp.nc$var[[i]]$name)
			}
		print("following varnames are given:")
		for (i in 1:length(varnames)){print(varnames[i])}	
		stop("you have to specify varname")
		}
		varname<-temp.nc$var[[1]]$name
	}
	
	# read name for lon and lat variabels from temp.nc
	if(is.null(lonname)){
		lonnames<-c("lon","longitude") # list of known lonnames
		lonname<-find.var(temp.nc,lonnames)[1]
	}
	if(is.null(latname)){
		latnames<-c("lat","latitude")	 # list of known latnames
		latname<-find.var(temp.nc,latnames)[1]
	}

	#Read out the data
	temp.time <- get.var.ncdf(temp.nc,"time")
	temp.data <-get.var.ncdf(temp.nc,varname)
	temp.lat <-get.var.ncdf(temp.nc,latname)
	temp.lon <-get.var.ncdf(temp.nc,lonname)
	
	#convert from missVal given values to NA
	temp.data[temp.data<=missVal[1]]<-NA
	temp.data[temp.data>=missVal[2]]<-NA


	##convert dates for yearly and monthly data
	# get informations about "time"-variable
	timevar<-as.numeric(find.var(temp.nc,"time")[2:3])
	unit.time<-temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$units
	diff.time<-max(diff(temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals))
	#diff.time<-temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals[[2]]-temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals[[1]]

	
	if(unit.time=="day as %Y%m%d.%f"){
		if(diff.time==100){
			year <- floor(temp.time/10000)
			temp.date <- year + (floor((temp.time-(year*10000))/100)-1)/12
		}else{
		if(diff.time==10000){
			temp.date<-temp.time%/%10000
		}else{
			if(min(diff(temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals))==1){
				d.year<-floor(temp.time/10000)
				
				reftime<-julday.own(floor(temp.time[1]/10000)*10000+101)
				d.day<-julday.own(temp.time)-reftime
				len<-length(temp.date)
				d.day[d.day>(len-1)]<-d.day[d.day>(len-1)]-len
				temp.date<-d.year+d.day/365
				
			}else{stop("time steps are not daily, monthly or yearly")}
		}}
	}else{
	if(unit.time=="hours since 1-1-1 00:00:0.0"|unit.time=="hours since 1-01-01 00:00"){
		if (diff.time==24){
			temp.date<-(chron(temp.time/24,origin=c(month=1,day=1,year=01)))
			d.year<-as.numeric(as.character(years(temp.date)))
			d.day<-as.numeric(temp.date-chron(paste("1/1/",years(temp.date),sep="")))
			temp.date<-d.year+d.day/365

		}else{
		temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))
		}
	}else{
	if(length(grep(glob2rx("days since ????-??-?? ??:??"),unit.time))){
		start.year<-as.numeric(sub("-..-.....:..","",sub("days since ","",unit.time)))
		start.mon<-as.numeric(sub("-.....:..","",sub("days since ....-","",unit.time)))
		start.day<-as.numeric(sub("...:..","",sub("days since ....-..-","",unit.time)))
		abs.start.day<-julday(start.mon,start.day,2001)-julday(1,1,2001)

		d.day<-(temp.time+abs.start.day)/365
		temp.date<-start.year+d.day
		
	}else{stop(paste("time format",unit.time,"not supported by read_data"))}
	}}


	#Sort the latitudes
	tmp<-sort(temp.lat,index.return=TRUE)
	temp.lat<-temp.lat[tmp$ix]
	temp.data<-temp.data[,tmp$ix,]
	
	#sort the longitudes
	temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360
	tmp<-sort(temp.lon,index.return=TRUE)
	temp.lon<-temp.lon[tmp$ix]
	temp.data<-temp.data[tmp$ix,,]


	return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}

find.var<-function(data.nc,searched_vars)
{
	for (i in 1:length(data.nc$var))
	{
		for (j in 1:length(data.nc$var[[i]]$dim))
		{
			if(is.element(data.nc$var[[i]]$dim[[j]]$name,searched_vars)){
				varname<-data.nc$var[[i]]$dim[[j]]$name
				return(c(varname,i,j))
			}
		}
	}	
}


julday.own<-function(x){
	year<-floor(x/10000)	
	month<-floor((x-year*10000)/100)
	day<-x-(year*10000+month*100)

	return(julday(month,day,year))
}


cor.pTs<-function(pTs,field,use="complete.obs",min.obs=30,debug=F)
{
        #bring both on the same time basis
        start<-max(start(pTs)[1],start(field)[1])
        end<-min(end(pTs)[1],end(field)[1])
        if (debug) print(paste("Common time period: ",start,end))
        pTs<-window(pTs,start,end)
        field<-window(field,start,end)

        if (class(field)[1]=="pField") #field correlation
        {
		    n.Time<-length(time(field))
		    class(field)<-"matrix"
		    index<-((n.Time-colSums(is.na(field)))>min.obs)
		    dat<-field[,index]
                result<-matrix(NA,1,ncol(field))
        
                tresult<-cor(dat,pTs,use=use)
                result[,index]<-tresult
                class(field)<-"pField"
                return(copyattr(result,field))
        }
        else return(cor(pTs,field,use=use))
}

cortest.pTs<-function(pTs,field,min.obs=30)
{
	#bring both on the same time basis
	start<-max(start(pTs)[1],start(field)[1])
	end<-min(end(pTs)[1],end(field)[1])
	print(paste("Common time period: ",start,end))
	pTs<-window(pTs,start,end)
	field<-window(field,start,end)
	

	#Filter out data which contain not enough timesteps
	n.Time<-length(time(field))
	class(field)<-"matrix"
      index<-((n.Time-colSums(is.na(field)))>min.obs)
      dat<-field[,index]
	result<-matrix(NA,2,ncol(field))

	tresult<-apply(dat,2,mycor.test,c(pTs))
	result[,index]<-tresult

	class(field)<-"pField"
	return(copyattr(result,field))
}


plotmap.pField<-function(plotdata,main=NULL,zlim=range(plotdata,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=NULL,FUN=NULL,shift=F,long=F,xlim=NULL,stype=0,sSub="",set.bg=NULL, ...)
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

        tmp<-plot.preparation(plotdata,shift,long)
        if (is.null(xlim))  xlim = range(tmp$lon, finite = TRUE)
        
        if (stype == 1) 
        {
filled.contour.own(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,xlim=xlim,color=palette,set.bg=set.bg,plot.title={
        title(main=main,sub=sSub); 
        addland(col="black");
        if (!is.null(FUN)) FUN(tmp$lon,tmp$lat,tmp$data)
       
        },plot.axes=axes.stype(gridcolor,tmp$lat,tmp$lon))
      } else
        filled.contour.own(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,xlim=xlim,color=palette,set.bg=set.bg,plot.title={
        title(main=main,sub=sSub);
        addland(col="black"); grid()
        if (!is.null(FUN)) FUN(tmp$lon,tmp$lat,tmp$data)
        
       }, ...)
}

filled.contour.own<-function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)),
    z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE),
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) -
        1), plot.title, plot.axes, key.title, key.axes, asp = NA,
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, set.bg = NULL,
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

    if(!is.null(set.bg)){
        usr<-par('usr')
        rect(usr[1],usr[3],usr[2],usr[4],col=set.bg)
    }

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


plot.Polygon<-function(sigline)
{
   col="grey60"
   for (i in 1:length(sigline)) polygon(sigline[[i]]$x,sigline[[i]]$y,angle=-30,density=30,col=col, border = col)
}

plot.sig<-function(plotmap,sigmap,levels=0.95,...)
{

	temp<-plot.preparation(sigmap)

	diff.lon<-max(diff(temp$lon))
	diff.lat<-max(diff(temp$lat))
	len.lon<-length(temp$lon)
	len.lat<-length(temp$lat)
	
	new.lon<-c(temp$lon[1]-diff.lon,temp$lon,temp$lon[len.lon]+diff.lon)
	new.lat<-c(temp$lat[1]-diff.lat,temp$lat,temp$lat[len.lat]+diff.lat)
	
	empty.row<-matrix(ncol=len.lat)
	empty.col<-matrix(nrow=len.lon+2)
	
	new.data<-rbind(empty.row,temp$data)
	new.data<-rbind(new.data,empty.row)
	new.data<-cbind(empty.col,new.data)
	new.data<-cbind(new.data,empty.col)
	
	new.data[is.na(new.data)]<-1
	siglines<-contourLines(new.lon,new.lat,1-new.data,levels=levels)
	plot(plotmap,FUN=plot.Polygon(siglines),...)
}


filter.pField<-function(field,Filter,f.time,...)
{
	result<-apply(field,2,"filter",Filter,...)
	newTime<-window(time(field),start(field)[1]+f.time,end(field)[1]-f.time)
	newField<-pField(NULL,newTime,getlat(field),getlon(field),paste(getname(field)),gethistory(field),date=FALSE)

	for (i in 1:(length(time(field))-2*f.time))
	{
		newField[i,]<-result[i+f.time,]
	}
	return(newField)
      #return(pField(result,time(data),getlat(data),getlon(data),paste(getname(data),"filtered"),gethistory(data)))
	# obiges geht nicht, weil pField nicht mit NA-Werten klarkommt...
}


