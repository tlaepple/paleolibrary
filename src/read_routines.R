### The following functions all belong to read_data...
### unit.X convert the timesteps dependent on the units given
### They are called from getTimeAxes, which is called by the main function read_data

### This is the first approach of splitting and cleaning up the "univeral reading function"
### There is still lots of work to be done

unit.1<-function(diff.time,temp.time,temp.nc,timevar,unit.time)    #day as %Y%m%d.%f
{
    if(diff.time<10000){
            year <- floor(temp.time/10000)
            temp.date <- year + (floor((temp.time-(year*10000))/100)-1)/12
          }else
{  #1000 or multiples are interpreted as years
            if((min(diff(temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals)) == diff.time))
            {
              temp.date<-temp.time%/%10000
            } else{
              if (min(diff(temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals))==1) {
                d.year<-floor(temp.time/10000)
                reftime<-julday.own(floor(temp.time[1]/10000)*10000+101)
                d.day<-julday.own(temp.time)-reftime
                len<-length(temp.date)
                d.day[d.day>(len-1)]<-d.day[d.day>(len-1)]-len
                temp.date<-d.year+d.day/365
              }else{stop("time steps are not daily, monthly or yearly")}
            }
 }
    return(temp.date)
}

unit.2<-function(diff.time,temp.time,temp.nc,timevar,unit.time)### ="hours since 1-1-1 00:00:0.0"|unit.time=="hours since 1-01-01 00:00
{
  if (diff.time==24){
              temp.date<-(chron(temp.time/24,origin=c(month=1,day=1,year=01)))
              d.year<-as.numeric(as.character(years(temp.date)))
              d.day<-as.numeric(temp.date-chron(paste("1/1/",years(temp.date),sep="")))
              temp.date<-d.year+d.day/365

            } else{
              temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))
            }
  return(temp.date)
}

unit.3<-function(diff.time,temp.time,temp.nc,timevar,unit.time) # days since ????-??-?? ??:??
{
    start.year<-as.numeric(sub("-..-.....:..","",sub("days since ","",unit.time)))
    start.mon<-as.numeric(sub("-.....:..","",sub("days since ....-","",unit.time)))
    start.day<-as.numeric(sub("...:..","",sub("days since ....-..-","",unit.time)))
    abs.start.day<-julday(start.mon,start.day,2001)-julday(1,1,2001)

    d.day<-(temp.time+abs.start.day)/365
    temp.date<-start.year+d.day
    return(temp.date)
}

unit.4<-function(diff.time,temp.time,temp.nc,timevar,unit.time) #    #days since ????-??-?? ??:??:??
{
    start.year<-as.numeric(sub("-..-.....:..:..","",sub("days since ","",unit.time)))
    start.mon<-as.numeric(sub("-.....:..:..","",sub("days since ....-","",unit.time)))
    start.day<-as.numeric(sub("...:..:..","",sub("days since ....-..-","",unit.time)))
    abs.start.day<-julday(start.mon,start.day,2001)-julday(1,1,2001)
    d.day<-(temp.time+abs.start.day)/365
    temp.date<-start.year+d.day
    return(temp.date)
}


unit.5<-function(diff.time,temp.time,temp.nc,timevar,unit.time) #  days since ???-??-?? ??:??:??
{

    start.year<-as.numeric(sub("-..-.....:..:..","",sub("days since ","",unit.time)))
    start.mon<-as.numeric(sub("-.....:..:..","",sub("days since ...-","",unit.time)))
    start.day<-as.numeric(sub("...:..:..","",sub("days since ...-..-","",unit.time)))

    temp.date <- as.vector(as.yearmon(chron(temp.time,origin=c(month=start.mon,day=start.day,year=start.year))))
    # cut after comma
    temp.date<-floor(temp.date)

    return(temp.date)
}






getTimeAxes<-function(temp.time,temp.nc)
#Input:
#temp.time: time vector in the orginal format
#temp.nc  : link to the netcdf file

#Output: Date vector in year.year
{


    timevar<-as.numeric(find.var(temp.nc,"time")[2:3])
    unit.time<-temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$units
#UNCLEAR; can we just use temp.time here?
    tval<-temp.nc$var[[timevar[1]]]$dim[[timevar[2]]]$vals
    if (length(tval)>1) diff.time<-max(diff(tval)) else diff.time=0


     if(unit.time=="day as %Y%m%d.%f") temp.date<-unit.1(diff.time,temp.time,temp.nc,timevar,unit.time) else
     if(unit.time=="hours since 1-1-1 00:00:0.0"|unit.time=="hours since 1-01-01 00:00")  temp.date<-unit.2(diff.time,temp.time,temp.nc,timevar,unit.time) else
     if(length(grep(glob2rx("days since ????-??-?? ??:??"),unit.time)))  temp.date<-unit.3(diff.time,temp.time,temp.nc,timevar,unit.time) else
     if (length(grep(glob2rx("days since ????-??-?? ??:??:??"),unit.time))) temp.date<-unit.4(diff.time,temp.time,temp.nc,timevar,unit.time) else
     if(length(grep(glob2rx("days since ???-??-?? ??:??:??"),unit.time))) temp.date<-unit.5(diff.time,temp.time,temp.nc,timevar,unit.time) else
               stop(paste("time format",unit.time,"not supported by read_data"))

     return(temp.date)
}



read_data<-function(FILENAME="",varname=NULL,name="",lonname=NULL,latname=NULL,timename=NULL,missVal=c(-1e+20,1e+20),bClim=FALSE,level=NULL,roundMonth=FALSE)
#Universal reading function for netcdf files; originally from Moritz Krieger
#Input:
#    FILENAME
#    varname
#optionaly: lonname,latname and timename, if these are not the standard names
#missVal: Values outside this range are set to NA
#bClim (TRUE): than it expects 12 values and gives them the months 1...12
#roundMonth (TRUE): rounds the dates to the next month; this is sometimes needed to create
#equidistant timesteps (as the real months have different lengths)

#level: level index; if a multilevel variable is choosen and level is given, this level is returned

#Output
#pField object containing the variable + lat,lon,time and the name
#26.10.11: The reading process was splitted into functions to keep the code more readable; still untested
{

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
        lonnames<-c("lon","longitude","lons") # list of known lonnames
        lonname<-find.var(temp.nc,lonnames)[1]
    }
    if(is.null(latname)){
        latnames<-c("lat","latitude","lats")     # list of known latnames
        latname<-find.var(temp.nc,latnames)[1]
    }

       if(is.null(timename)){
        timenames<-c("time","t")     # list of known timenames
        timename<-find.var(temp.nc,timenames)[1]
    }


    #Read out the data
    temp.time <- get.var.ncdf(temp.nc,timename)
    temp.data <-get.var.ncdf(temp.nc,varname)
    temp.lat <-get.var.ncdf(temp.nc,latname)
    temp.lon <-get.var.ncdf(temp.nc,lonname)
       #convert from missVal given values to NA
    temp.data[temp.data<=missVal[1]]<-NA
    temp.data[temp.data>=missVal[2]]<-NA

    if (!is.null(level)) #User provided a level
    {
        #if only one timestep; than level = third dimension,
        if (length(temp.time)<2)  temp.data<-temp.data[,,level] else temp.data<-temp.data[,,level,]
        warning("Be carful; the level reading might not yet always work")
    }


    if (bClim & (length(temp.time))!=12) stop("bClim was choosen but != 12 timesteps in the dataset")

    if (!bClim) temp.date<-getTimeAxes(temp.time,temp.nc) else temp.date=1:12


    if (length(dim(temp.data))==3)
    {

    #Sort the latitudes
    tmp<-sort(temp.lat,index.return=TRUE)
    temp.lat<-temp.lat[tmp$ix]
    temp.data<-temp.data[,tmp$ix,]

       #sort the longitudes
    temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360
    tmp<-sort(temp.lon,index.return=TRUE)
    temp.lon<-temp.lon[tmp$ix]
    temp.data<-temp.data[tmp$ix,,]
  }  else #temp.time has only one element
      {

    if (length(dim(temp.data))!=2) stop("three dimensions were found but only one timestep; if the file includes levels, please provide the level to read")
          #Sort the latitudes

    tmp<-sort(temp.lat,index.return=TRUE)
    temp.lat<-temp.lat[tmp$ix]
    temp.data<-temp.data[,tmp$ix]

       #sort the longitudes
    temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360
    tmp<-sort(temp.lon,index.return=TRUE)
    temp.lon<-temp.lon[tmp$ix]
    temp.data<-temp.data[tmp$ix,]
      }


if (roundMonth)
{
    warning("each timestep is rounded to the next full month")
    temp.date<-round(temp.date*12)/12
}


#Check if multiple timesteps are contained; For single fields (e.g. Trends) return the single field with time 1
   if (length(temp.date)>1) return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME)) else
    return(pField(temp.data,1,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
 }

read.rodgers.nao<-function(FILENAME="/home/tlaepple/data/optimalforecast/Rogers_NAO_monthly_normalized_anomalies.txt")
{
temp<-read.table(FILENAME,skip=3,na.strings="-9999.0")
nao.data<-t(as.matrix(temp[,2:13]))
all<-c(nao.data)
date<-rep(temp[,1],each=12)+(0:11)/12+0.01
nao<-pTs(all,date,0,0,name="NAO Rogers")
return(nao)
}

read.angsmalik<-function(FILENAME="c:/data/coupled/Angsmallik.txt")
{
temp<-read.csv(FILENAME,skip=1,na.strings="99.9",fill=T,header=F)
data<-t(as.matrix(temp[,1:12]))
data[data==99.9]<-NA
all<-c(data)
date<-rep(1895:2003,each=12)+(0:11)/12
data<-pTs(all,date,0,0,name="NAO Rogers")
return(data)
}

read_solarforcing<-function(FILENAME="~/data/sonne/forcing/solar.dat")
{
data<-read.table(FILENAME,header=TRUE)
data<-pTs(data[,2],data[,1],NA,NA,name="SF: W/m2",FILENAME)
return(data)
}

read_gasforcing<-function(FILENAME="~/data/sonne/forcing/gases.dat")
{
data<-read.table(FILENAME,header=TRUE)
data<-pTs(data[,2],data[,1],NA,NA,name="CO2: ppm",FILENAME)
return(data)
}

#Routines for reading fields
#Read routines must supply pField types with lats and lons in a strictly increasing order..




#uses ncdf,chron
read_ncep <- function(FILENAME="",varname="",name="")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
#temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))
temp.date<-temp.time%/%10000
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}



read_mann<-function(FILENAME="c:/data/NAO_SOI/mannSOI.dat",index=2)
{
data<-read.table(FILENAME,header=TRUE,na.strings="-99.99")
data<-pTs(data[,index],data[,1],NA,NA,name="SOI",FILENAME)
return(data)
}


read_precipdata<-function(FILENAME="c:/data/norel/RR_LOCID000254_precip.txt",skip=20,lat=69.21,lon=360-51.1,name="prec Illuisat",na.pad=TRUE)
{
data<-read.csv(FILENAME,skip=skip,na.strings="-9999")

#data gymnastics
date<-data[,2]
d.year<-date%/%10000
d.month<-(date-d.year*10000)%/%100
d.day<-(date-d.year*10000-d.month*100)

date<-chron(paste(d.month,"/",d.day,"/",d.year))
d.day<-as.numeric(date-chron(paste("1/1/",years(date),sep="")))

index.without366<-(d.day<365)
temp.date<-d.year+d.day/365
temp.date<-temp.date[index.without366]

temp.data<-data[index.without366,3]
#Fill NA with zeros
if (na.pad) temp.data[is.na(temp.data)]<-0
result<-pTs(temp.data,temp.date,lat,lon,name,FILENAME)
return(result)
}


read_ncep_yr<- function(FILENAME,varname)
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
#temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))

temp.date<-temp.time%/%10000
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}


#uses ncdf,chron... skips the 366's day...-> every year = 365 days
read_ncep_day<- function(FILENAME,varname)
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
#temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))

temp.date<-(chron(temp.time/24,origin=c(month=1,day=1,year=01)))

d.year<-as.numeric(as.character(years(temp.date)))
d.day<-as.numeric(temp.date-chron(paste("1/1/",years(temp.date),sep="")))

index.without366<-(d.day<365)
temp.date<-d.year+d.day/365
temp.date<-temp.date[index.without366]

temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")


temp.data <- temp.data[,,index.without366]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}



#uses ncdf,skips the 366's day...-> assumes that the file contains day 1:365 or 366
read_ncep.clim.day<- function(FILENAME,varname)
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
temp.date<-1:365

temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")


temp.data <- temp.data[,,1:365]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}










#read monthly ncep dataset
#uses ncdf,chron
read_ipcc.mon <- function(FILENAME="C:/data/NCEP/DJFslp.nc",varname="slp",sy=1850)
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
temp.date <- as.vector(as.yearmon(chron(temp.time,origin=c(month=1,day=1,year=1850))))

temp.date<-sy+1/24+1/12*(1:length(temp.time))

temp.data <-get.var.ncdf(temp.nc,varname)

 temp.data[temp.data>1e5]<-NA

temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}


#read monthly ncep dataset
#uses ncdf,chron
read_ncep.mon <- function(FILENAME="C:/data/NCEP/DJFslp.nc",varname="slp")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}


#read sonne dataset (Martin Stendel's model output converted with cdo)
read_sonne <- function(FILENAME="~/data/sonne/DJF139.nc",GRIDNAME="~/data/t42.nc",varname="var139",name=NULL)
{

temp.nc = open.ncdf(FILENAME)

if (is.null(name)) name<-FILENAME

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
year <- floor(temp.time/10000)
temp.date <- year #+ (floor((temp.time-(year*10000))/100)-1)/12
grid.nc =  open.ncdf(GRIDNAME)
temp.data <-get.var.ncdf(temp.nc,varname)

temp.lat <-get.var.ncdf(grid.nc,"lat")
temp.lon <-get.var.ncdf(grid.nc,"lon")

tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
#transform into 2D first
dim(temp.data)<-c(length(temp.lon),length(temp.lat),length(temp.date))
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}

#read kaplan dataset
read_kaplan <- function(FILENAME="",varname="ssta")
{

temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
temp.data <-get.var.ncdf(temp.nc,varname)

year <- floor(temp.time/10000)
temp.date <- year #+ (floor((temp.time-(year*10000))/100)-1)/12


temp.lat <-get.var.ncdf(temp.nc,"Y")
temp.lon <-get.var.ncdf(temp.nc,"X")
temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]



return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name="Kaplan SST",history=FILENAME))
}


#read kaplan dataset
read_kaplan_monthly <- function(FILENAME="",varname="sst")
{

temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
n<-length(temp.time)

temp.data <-get.var.ncdf(temp.nc,varname)
temp.date <- 1856+1/24+(1:n)/12


temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")
temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]



return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name="Kaplan SST",history=FILENAME))
}




load_huascara<-function(FILENAME="c:/data/huascara/huascara.txt")
{
	data<-read.table(FILENAME,header=TRUE)
	n<-length(data[,1])
	result<-pTs(data[n:1,2:5],data[n:1,1],-999,-999,c("dO18","particle","NO3-","dO18_2"),FILENAME)
	return(result)

}




read.accum <- function(filename,name,lat,lon) #function to read accumulation data
  {
    data<-read.table(filename)
    return(pTs(data[,2],data[,1],lat,lon,name,filename))

 }

read.had.monthly<-function(FILENAME="",varname="temp")
{

temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
n<-length(temp.time)

temp.data <-get.var.ncdf(temp.nc,varname)
temp.date <- 1870+(0:(n-1))/12



temp.lat <-get.var.ncdf(temp.nc,"latitude")
temp.lon <-get.var.ncdf(temp.nc,"longitude")
temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

temp.data[temp.data>1e20]<-NA


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name="HADI
SST",history=FILENAME))
}



read.had.annual<-function(FILENAME="",varname="temp",latname="latitude",lonname="longitude")
{

temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
temp.date<-floor(temp.time/10000)

temp.data <-get.var.ncdf(temp.nc,varname)


temp.lat <-get.var.ncdf(temp.nc,latname)
temp.lon <-get.var.ncdf(temp.nc,lonname)
temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

temp.data[temp.data>1e20]<-NA


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name="HADI
SST",history=FILENAME))
}




read.mon.ecmwf <- function(FILENAME="/home/tlaepple/data/optimalforecast/ecmwf/gph500.mon.nc",varname="z")
  {
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")

year<-floor(temp.time/10000)
month<-floor((temp.time-year*10000)/100)
temp.date<-year+(month-1)/12

temp.data <-get.var.ncdf(temp.nc,varname)


temp.lat <-get.var.ncdf(temp.nc,"latitude")
temp.lon <-get.var.ncdf(temp.nc,"longitude")
temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]

temp.data[temp.data>1e20]<-NA


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name="HADI
SST",history=FILENAME))
  }





read.clim<-function(FILENAME="",varname="",lonname="lon")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.date<-1:12
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,lonname)


temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}



read.ice.had.annual<-function(FILENAME="",varname="")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.date<-1870:2007
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")


temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}





read.ice.had.sd<-function(FILENAME="",varname="")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.date<-1
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")


temp.lon[temp.lon<0]<-temp.lon[temp.lon<0]+360

#Sort longitudes
tmp<-sort(temp.lon,index.return=TRUE)
temp.lon<-temp.lon[tmp$ix]
temp.data<-temp.data[tmp$ix,]

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=FILENAME,history=FILENAME))
}



read.mld.kara <- function(FILENAME="C:/data/MLD/mld0p8t.woa.global.mc.1994.ncdf",varname="MIXED_LAYER_DEPTH",name="")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"Time")
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"Latitude")
temp.lon <-get.var.ncdf(temp.nc,"Longitude")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.time,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}





#uses ncdf,chron
read_ncep.clim <- function(FILENAME="C:/data/NCEP/DJFslp.nc",varname="slp",name="",lonname="lon",latname="lat")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.time <- get.var.ncdf(temp.nc,"time")
#temp.date <- as.vector(as.yearmon(chron(temp.time/24,origin=c(month=1,day=1,year=01))))
temp.date<-1:12
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,latname)
temp.lon <-get.var.ncdf(temp.nc,lonname)

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix,]


return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}





#uses ncdf,chron
read_one <- function(FILENAME="C:/data/NCEP/DJFslp.nc",varname="slp",name="")
{
temp.nc = open.ncdf(FILENAME)

#Read out the data
temp.data <-get.var.ncdf(temp.nc,varname)
temp.lat <-get.var.ncdf(temp.nc,"lat")
temp.lon <-get.var.ncdf(temp.nc,"lon")

#Sort the latitudes
tmp<-sort(temp.lat,index.return=TRUE)
temp.lat<-temp.lat[tmp$ix]
temp.data<-temp.data[,tmp$ix]


return(pField(temp.data,1,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}



read_er.clim<-function(FILENAME,varname="",name=""){
        temp.nc = open.ncdf(FILENAME)

        # read name for lon and lat variabels from temp.nc
        name.lon<-"lon_2"
        name.lat<-"lat"

        if(varname==""){
                varname<-temp.nc$var[[1]]$name
                if(temp.nc$nvars>1){warning("no varname given and more then one variable")}
        }

        #Read out the data
        temp.time <- get.var.ncdf(temp.nc,"time")
        temp.data <-get.var.ncdf(temp.nc,varname)
        temp.lat <-get.var.ncdf(temp.nc,name.lat)
        temp.lon <-get.var.ncdf(temp.nc,name.lon)



                        temp.date <- 1:12



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




read_er<-function(FILENAME,varname="",name=""){
        temp.nc = open.ncdf(FILENAME)

        # read name for lon and lat variabels from temp.nc
        name.lon<-"lon_2"
        name.lat<-"lat"

        if(varname==""){
                varname<-temp.nc$var[[1]]$name
                if(temp.nc$nvars>1){warning("no varname given and more then one variable")}
        }

        #Read out the data
        temp.time <- get.var.ncdf(temp.nc,"time")
        temp.data <-get.var.ncdf(temp.nc,varname)
        temp.lat <-get.var.ncdf(temp.nc,name.lat)
        temp.lon <-get.var.ncdf(temp.nc,name.lon)



                        year <- floor(temp.time/10000)
                        temp.date <- year



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

