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


getx<-function(x,i) return(x[i]) #Trivial function to retrun the ith element of a vector x

read_data<-function(FILENAME="",varname=NULL,name="",lonname=NULL,latname=NULL,timename=NULL,
                    bIndexTime=FALSE,levelname=NULL,missVal=c(-1e+20,1e+20),bClim=FALSE,level=NULL,roundMonth=FALSE)
#Universal reading function for netcdf files; originally from Moritz Krieger
#Input:
#    FILENAME
#    varname
#optionaly: lonname,latname and timename, if these are not the standard names
#missVal: Values outside this range are set to NA
#bIndexTime (TRUE): ignores the time variable and uses the sequential index as time
#roundMonth (TRUE): rounds the dates to the next month; this is sometimes needed to create
#equidistant timesteps (as the real months have different lengths)
#levelname: must be given if a level shopuld be extracted
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
      # todo: analoque for levels

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
        level.index <- find.dimposition(temp.nc,varname,levelname)
        if(is.null(level.index)) stop('Level name does not exist.')
        dim.not.level<-(1:length(dim(temp.data)))[-1*level.index] #Return all dimension indices except level.index
        temp.data<-apply(temp.data,dim.not.level,getx,i=level) #Return the array from selected level 
         warning("Be carful; the level reading might not yet always work")
    }


    if (!bIndexTime) temp.date<-getTimeAxes(temp.time,temp.nc) else temp.date=seq(temp.time)


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
