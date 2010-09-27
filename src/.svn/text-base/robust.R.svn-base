### Robust statistics; 22.6.
# Aim: Robust mean and variability+ Outlier detection
# Date:
# Documentation: Mudelsee 2006
#Short note
#CLIM-X-DETECT: A Fortran 90 program for robust
#detection of extremes against a time-dependent background in
#climate reconstruction


## running Median and MAD 
medsmooth<-function(x,m)
{
#Median von i0.i.i1
N<-length(x)
rmedian<-rmad<-rmean<-vector()
	for (i in 1:N)
	{
		i0<-max(i-m,1)
		i1<-min(i+m,N)
		part<-x[i0:i1]
		rmedian[i]<-median(part)
		rmad[i]<-mad(part)
		rmean[i]<-mean(part)
	}


result<-list(rmedian=rmedian,rmad=rmad,rmean=rmean,m=m)
class(result)<-"msmooth"
return(result)
}

#Testcode for medsmooth
#test<-rnorm(100)
#test[30]<-20
#test[66:67]<-8
#test[80]<--4
#plot(test)
#mpr<-medsmooth(test,5)
#lines(mpr$rmedian,col="red")
#lines(mpr$rmedian+mpr$rmad*3.5,col="red",lty=2)
#lines(mpr$rmad,col="green")
#lines(mpr$rmean,col="blue")
