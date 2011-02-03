
# load all needed packages
library(chron) #includes some time handling functions
library(clim.pact)  #includes the coastline drawing function
library(zoo)   #some timeseries function
library(ncdf)  #netcdf handling
library(lattice) #wireframe plot
library(Rwave) #Wavelet analysis
#library(sowas)
                                        #library(spatstat) #Spatial statistics


source(paste(path,'basis.R',sep=""))
source(paste(path,'proxytype.R',sep=""))        #type definition
source(paste(path,'selspace.pField.R',sep=""))  #selection of areas and points
source(paste(path,'mat.pField.R',sep=""))       #EOF's, correlation etc...

source(paste(path,'plotlib.R',sep=""))       #plotting routines
source(paste(path,'read_routines.R',sep="")) #reading routines
source(paste(path,'indices.R',sep=""))       #standard climate indices
source(paste(path,'roll.pField.R',sep=""))   #running/rolling functions

source(paste(path,'noise.pField.R',sep=""))  #surrogate timeseries for MonteCarlo

source(paste(path,'filter.R',sep=""))        #lowpass and bandpass filter caclulation
source(paste(path,'season.R',sep=""))                #season averages/
source(paste(path,'insolation.R',sep=""))
source(paste(path,'zonalmean.R',sep=""))                              #Insolation Code /
source(paste(path,'coherency.R',sep=""))
source(paste(path,'confspec.R',sep=""))

source(paste(path,'paleo.symbols.R',sep=""))



source(paste(path,'ownfunctions.R',sep=""))                #put your own functions /
							                  #changed functions here


source(paste(path,'selspace.interp.R',sep=""))


source(paste(path,'robust.R',sep=""))

source(paste(path,'confspec.R',sep=""))  #Includes spectral estimation + significance testing and conficence intervals

source(paste(path,'shortcuts.R',sep=""))  #Includes shortcuts as first and last to access the first and last element of a vector


