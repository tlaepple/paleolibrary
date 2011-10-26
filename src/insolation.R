
Tyear<-365.2422 #length of tropical year in days

#The following calculations were based on the book
#Grossman; The sheer joy of celestial mechanics
#Orbits under the inverse square law

calcE<-function(M,e,tol=1e-8)
#Calculate the Eccentricity anomaly E from the mean anomaly M, iterativly
#inversion of Keplers formula: E-e*sin(E)=M  (188)
{
  E<-vector()
  for (i in 1:length(M))
  {
    Etemp<-M[i]
    ratio<-1
    while (abs(ratio)>tol)
    {
        f.E <- Etemp - e*sin(Etemp) -M[i]
        f.Eprime <- 1-e*cos(Etemp)
        ratio <- f.E/f.Eprime
        if (abs(ratio) > tol) Etemp<-Etemp-ratio else E[i]<-Etemp

    }
}
    return(E)

}



calcf<-function(E,e)
{

#Calculate the true anomaly f (theta) from the eccentricity anomaly
 #true anomaly = angle between the point and perihelion
#book: (181)
    temp<-atan(sqrt((1+e)/(1-e))*tan(0.5*E))*2
temp[temp<0]<-temp[temp<0]+2*pi
    return(temp)
}


t.from.f<-function(ecc,f)
#input : f = true angle from perihelion
#output: t = days from perihelion, calculated using the mean angle M
{
    #Calculate the eccentricity anomaly E
     E<-atan(tan(0.5*f)/sqrt((1+ecc)/(1-ecc)))*2 #(181)
    #Keplers equation to get the mean anomaly (188)
     M <-  E - ecc*sin(E)

    #Convert into days
    n = 2*pi/Tyear
    t <- M/n #182
    if (t<0) t<-t+Tyear
    return(t)
}


f.from.t<-function(ecc,t)
{
#input: t = days from perihelion
#output : f = true angle from perihelion (via the mean angle M)
      n = 2*pi/Tyear
    M<-t*n
    E<-calcE(M,ecc)
    f<-calcf(E,ecc)
    return(f)
}






daily_insolation_param<-function(lat,day,ecc,obliquity,long_perh,day_type=1,alpha=0,T.alpha=80)
{
### Insolation, converted and adapted from the Eisenmann and Huybers,2006
#Matlab Code, based on Berger 1991
# Usage:
#   Fsw = daily_insolation(lat,day,eccentricity,obliquity,long_perh)
#
# Optional inputs/outputs:
#   [Fsw, ecc, obliquity, long_perh] = daily_insolation(kyear,lat,day,day_type)
#
# Description:
#   Computes daily average insolation as a function of day and latitude at
#   any point during the past 5 million years.
#
# Inputs:
#   kyear:    Thousands of years before present (0 to 5000).
#   lat:      Latitude in degrees (-90 to 90).
#   day:      Indicator of time of year; calendar day by default.
#   day_type: Convention for specifying time of year (+/- 1,2) [optional].
#     day_type=1 (default): day input is calendar day (1-365.24), where day 1
#       is January first.  The calendar is referenced to the vernal equinox
#       which always occurs at day 80.
#     day_type=2: day input is solar longitude (0-360 degrees). Solar
#       longitude is the angle of the Earth's orbit measured from spring
#       equinox (21 March). Note that calendar days and solar longitude are
#       not linearly related because, by Kepler's Second Law, Earth's
#       angular velocity varies according to its distance from the sun.
#    day_type=3; Custom calendar
#    that alpha (angle since equinox) is reached at day T.alpha
# this uses an iterative solution of Keplers equation and is therefore
# much slower than the day_type=1 which uses a trigeometric approximation
#of Kepler's equation

# Output:
#   Fsw = Daily average solar radiation in W/m^2.
#   omega = Angle in degree since the spring equinox
#
# This script contains orbital parameter data for the past 50000 years
#   from Berger and Loutre (1991).
#
# Detailed description of calculation:
#   Values for eccentricity, obliquity, and longitude of perihelion for the
#   past 5 Myr are taken from Berger and Loutre 1991 (data from
#   ncdc.noaa.gov). If using calendar days, solar longitude is found using an
#   approximate solution to the differential equation representing conservation
#   of angular momentum (Kepler's Second Law).  Given the orbital parameters
#   and solar longitude, daily average insolation is calculated exactly
#   following Berger 1978.
#
# References:
#   Berger A. and Loutre M.F. (1991). Insolation values for the climate of
#     the last 10 million years. Quaternary Science Reviews, 10(4), 297-317.
#   Berger A. (1978). Long-term variations of daily insolation and
#     Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
#     2362-2367.
#   Grossmann (1937). The sheer joy of celestial mechanics
#
# Authors:
#   Ian Eisenman and Peter Huybers, Harvard University, August 2006
#   eisenman@fas.harvard.edu
#   This file is available online at
#   http://deas.harvard.edu/~eisenman/downloads
#   Translated into R by Thomas Laepple
# Suggested citation:
#   P. Huybers and I. Eisenman, 2006. Integrated summer insolation
#   calculations. NOAA/NCDC Paleoclimatology Program Data
#   Contribution #2006-079.
#

# === Get orbital parameters ===


    epsilon<-obliquity*pi/180
	omega<-long_perh*pi/180

# === Calculate insolation ===
lat=lat*pi/180 # latitude

# lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
if (day_type ==1) # calendar days
{
    # estimate lambda from calendar day using an approximation from Berger 1978 section 3
    delta_lambda_m <- (day-80)*2*pi/365.2422;  ## das wäre lambda bei gleich langen Tagen

    beta<-(1-ecc^2)^(1/2)

    lambda_m0<-(-2)*( (1/2*ecc+1/8*ecc^3)*(1+beta)*sin(-omega)-1/4*ecc^2*(1/2+beta)*sin(-2*omega)+1/8*ecc^3*(1/3+beta)*(sin(-3*omega)))

    lambda_m<-lambda_m0+delta_lambda_m

    lambda<-lambda_m+(2*ecc-1/4*ecc^3)*sin(lambda_m-omega)+(5/4)*ecc^2*sin(2*(lambda_m-omega))+(13/12)*ecc^3*sin(3*(lambda_m-omega))


} else if (day_type==2) #solar longitude (1-360)
    {lambda=day*2*pi/360 # lambda=0 for spring equinox
} else if (day_type==3)
{
    #t.from.f(eccc,alpha-omega) = days since Perihelion in which
    #the angle alpha, measured from the equinox is reached

    #tau
    tau = -1*t.from.f(ecc,alpha-omega) + T.alpha
    #angle since equinox given the days with the calendar fixed using
    #the condition that alpha (angle since equinox) is reached at day T.alpha
    lambda<-f.from.t(ecc,day-tau)+omega



}else stop('Error: invalid day_type')


So<-1365 # solar constant (W/m^2)
delta<-asin(sin(epsilon)*sin(lambda)) # declination of the sun
Ho<-acos(-tan(lat)*tan(delta)) # hour angle at sunrise/sunset

# no sunrise or no sunset: Berger 1978 eqn (8),(9)

Ho[(abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta > 0 )]  <- pi
Ho[( abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta <= 0 )] <- 0

# Insolation: Berger 1978 eq (10)
Fsw=So/pi*(1+ecc*cos(lambda-omega))^2 /(1-ecc^2)^2 * ( Ho*sin(lat)*sin(delta) + cos(lat)*cos(delta)*sin(Ho))

return(list(Fsw=Fsw,ecc=ecc,obliquity=obliquity,long_perh=long_perh,lambda=lambda/2/pi*360 ))
}



annual_insolation<-function(kyear,lat)
{
    result<-vector()
    for (i in 1:length(kyear)) result[i]<-mean(daily_insolation(kyear[i],lat,1:365)$Fsw)
    return(result)
}


# This script contains orbital parameter data for the past 50000 years
#   from Berger and Loutre (1991).
#
# Detailed description of calculation:
#   Values for eccentricity, obliquity, and longitude of perihelion for the
#   past 5 Myr are taken from Berger and Loutre 1991 (data from
#   ncdc.noaa.gov). If using calendar days, solar longitude is found using an
#   approximate solution to the differential equation representing conservation
#   of angular momentum (Kepler's Second Law) (day_type=1), or an iterative
#   solution of Keplers Equation (day_type=3)
#
#    Given the orbital parameters
#   and solar longitude, daily average insolation is calculated exactly
#   following Berger 1978.
#
# References:
#   Berger A. and Loutre M.F. (1991). Insolation values for the climate of
#     the last 10 million years. Quaternary Science Reviews, 10(4), 297-317.
#   Berger A. (1978). Long-term variations of daily insolation and
#     Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
#     2362-2367.

daily_insolation<-function(kyear,lat,day,day_type=1,fast=T,T.alpha=80,alpha=0)
{
### Computes the daily average insolation, converted and adapted from the Eisenmann and Huybers,2006
#Matlab Code
# Orbital parameters from Berger and Loutre 1991
# calculation of daily average insolation following Berger A. (1978).
# Iterative solution for the Kepler equation, following Grossman; "The sheer joy of celestial mechanics,
# Orbits under the inverse square law"
# Inputs:
#   kyear:    Thousands of years before present (0 to 5000).
#   lat:      Latitude in degrees (-90 to 90).
#   day:      Indicator of time of year; calendar day by default.
#   day_type: Convention for specifying time of year (1,2,3) [optional].
#     day_type=1 (default): day input is calendar day (1-365.24), where day 1
#       is January first.  The calendar is referenced to the vernal equinox
#       which always occurs at day 80.
#     day_type=2: day input is solar longitude (0-360 degrees). Solar
#       longitude is the angle of the Earth's orbit measured from spring
#       equinox (21 March). Note that calendar days and solar longitude are
#       not linearly related because, by Kepler's Second Law, Earth's
#       angular velocity varies according to its distance from the sun.
#  day_type=3; Custom calendar
#    day input is calendar day (1-365.24), where day 1
#    is January first.  The calendar is referenced to the solar longitude alpha
#    (0..2pi, angle of the earth orbit from spring equinox)
#       which always occurs at day T.alpha.
# this uses an iterative solution of Keplers equation and is therefore
# much slower than the day_type=1 which uses a trigeometric approximation
#of Kepler's equation

# Output:
#   Fsw = Daily average solar radiation in W/m^2.
#   omega = Angle in degree since the spring equinox
#
# Authors of original MATLAB version :
#   Ian Eisenman and Peter Huybers, Harvard University, August 2006
#   eisenman@fas.harvard.edu
#   This file is available online at
#   http://deas.harvard.edu/~eisenman/downloads
#
#Translation into R and extension to arbitray calendars; Thomas Laepple, 2008
#

# === Get orbital parameters ===

   if (fast) temp<-orbital_parameters_fast(kyear) else  temp<-orbital_parameters(kyear)
    ecc<-temp$ecc
    epsilon<-temp$epsilon
    omega<-temp$omega

# For output of orbital parameters
obliquity<-epsilon*180/pi
long_perh<-omega*180/pi

# === Calculate insolation ===
lat=lat*pi/180 # latitude

# lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
if (day_type ==1) # calendar days
{
    # estimate lambda from calendar day using an approximation from Berger 1978 section 3
    delta_lambda_m <- (day-80)*2*pi/365.2422;
    beta<-(1-ecc^2)^(1/2)
    lambda_m0<-(-2)*( (1/2*ecc+1/8*ecc^3)*(1+beta)*sin(-omega)-1/4*ecc^2*(1/2+beta)*sin(-2*omega)+1/8*ecc^3*(1/3+beta)*(sin(-3*omega)))
    lambda_m<-lambda_m0+delta_lambda_m
    lambda<-lambda_m+(2*ecc-1/4*ecc^3)*sin(lambda_m-omega)+(5/4)*ecc^2*sin(2*(lambda_m-omega))+(13/12)*ecc^3*sin(3*(lambda_m-omega))
} else if (day_type==2) #solar longitude (1-360)
    {lambda=day*2*pi/360 # lambda=0 for spring equinox
} else if (day_type==3)
{
    #t.from.f(eccc,alpha-omega) = days since perihelion in which
    #the angle alpha, measured from the equinox is reached

    #tau
    tau = -1*t.from.f(ecc,alpha-omega) + T.alpha

    #angle since equinox given the days with the calendar fixed using
    #the condition that alpha (angle since equinox) is reached at day T.alpha
    lambda<-f.from.t(ecc,day-tau)+omega



}else stop('Error: invalid day_type')

So<-1365 # solar constant (W/m^2)
delta<-asin(sin(epsilon)*sin(lambda)) # declination of the sun
Ho<-acos(-tan(lat)*tan(delta)) # hour angle at sunrise/sunset

# no sunrise or no sunset: Berger 1978 eqn (8),(9)

Ho[(abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta > 0 )]  <- pi
Ho[( abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta <= 0 )] <- 0

# Insolation: Berger 1978 eq (10)
Fsw=So/pi*(1+ecc*cos(lambda-omega))^2 /(1-ecc^2)^2 * ( Ho*sin(lat)*sin(delta) + cos(lat)*cos(delta)*sin(Ho))

return(list(Fsw=Fsw,ecc=ecc,obliquity=obliquity,long_perh=long_perh,lambda=lambda/2/pi*360 ))
}




orbital_parameters_fast<-function(kyear)
{
if ( length(orbital_global$ecc) != 50001) halt("orbital parameters not initalized")

return(list(ecc=orbital_global$ecc[kyear*10+1],epsilon=orbital_global$epsilon[kyear*10+1],omega=orbital_global$omega[kyear*10+1]))

}

##############
orbital_parameters<-function(kyear,FILEDATA="e:/data/insolation/data/ins_data.txt")
{
# === Load orbital parameters (given each kyr for 0-5Mya) ===
# Load the matrix contains data from Berger and Loutre (1991),
#   downloaded as ORBIT91 from ncdc.noaa.gov

m<-read.table(FILEDATA)

kyear0<- -1*m[,1] # kyears before present for data (kyear0>=0);
ecc0<-m[,2] # eccentricity

# add 180 degrees to omega (see lambda definition, Berger 1978 Appendix)
omega0<-m[,3]+180; # longitude of perihelion (precession angle)
omega0=unwrap(omega0*pi/180)*180/pi # remove discontinuities (360 degree jumps)

epsilon0=m[,4] # obliquity angle

# Interpolate to requested dates
ecc=spline(x=kyear0,y=ecc0,xout=kyear)$y
omega=spline(x=kyear0,y=omega0,xout=kyear)$y * pi/180
epsilon=spline(x=kyear0,y=epsilon0,xout=kyear)$y * pi / 180

return(list(ecc=ecc,epsilon=epsilon,omega=omega))
}

#Q = unwrap(P) corrects the radian phase angles in array P by adding multiples of ±2pi when absolute jumps between consecutive array elements are greater than pi radians.
# based on http://ccrma.stanford.edu/~jos/sasp/Matlab_listing_unwrap_m.html
unwrap<-function(p)
{

  N = length(p)
  up = rep(0,N)
  pm1 = p[1]
  up[1] = pm1
  po = 0
  thr = pi
  pi2 = 2*pi
  for (i in 2:N)
  {
    cp <- p[i] + po
    dp <- cp-pm1
    pm1 <- cp
    if (dp>=thr)
    {
      while (dp>=thr)
	{
        po <- po - pi2
        dp <- dp - pi2
      }
   }
    if (dp<=((-1)*thr))
    {
      while (dp<=thr)
	{
        po <- po + pi2
        dp <- dp + pi2
      }
    }
    cp = p[i] + po
    pm1 = cp
    up[i] = cp
  }
return(up)
}


#Validation

#plot(june.65N)
#june.65N.new<-vector()
#for (i in 1:10000)
#{
#	june.65N.new[i]<-daily_insolation(kyear=i/10,lat=65,day=172)$Fsw
#}
#
#
#
#lines(-(1:10000/10),june.65N.new,col="red")



tlag<-function(data,ilag)
{
	temp<-rep(data,3)
	return(temp[366:730-ilag])
}



### Cover functions to return the daily insolation, either using a classical calendar (aligned with March21) or using an alignment with Dec21 summer solstice
ins.march21<-function(kyear,LAT)
{
	return(daily_insolation(kyear,LAT,1:365)$Fsw)
}

ins.dec21<-function(kyear,LAT)
{
	r<-daily_insolation(kyear,LAT,1:365)  #Eigentlich passt 356 besser
	shift<-355-which.min(abs(r$lambda-270))	 #270    #dann entprechend hinschieben
	r<-tlag(r$Fsw,shift)
	return(r)
}



ins.dec21.param<-function(ecc,obliquity,long_perh,LAT)
{
	r<-daily_insolation_param(LAT,1:365,ecc,obliquity,long_perh)  #Eigentlich passt 356 besser
	shift<-355-which.min(abs(r$lambda-270))	     #dann entprechend hinschieben1
	r<-tlag(r$Fsw,shift)
	return(r)
}



orbital_global<-orbital_parameters((0:50000)/10,paste(path,"ins_data.txt",sep=""))



