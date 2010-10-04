# Aim:
# Date:
# Documentation:

#Changes:


### Über die ChiSq Verteilung und die Freiheitsgrade, kann die Varianz berechnet werden
## Das mittlere Spektrum kann über ein AR1 Prozess, der nach einem Median smoother angewandt wird, gerschätzt werden


specred<-function(S0,a0,omega,dt)
{
	fN=0.5/dt #Nyquist frequency
	return((S0*(1-a0^2)/(1-2*a0*cos(pi*omega/fN)+a0^2))*dt)
}



#sd value of AR1 power spectrum (Eq 5), given the sd of the noise and a0
getS0<-function(sigma,a0) return(sigma^2/(1-a0^2))


#Fitting a AR1 process in the spectral domain, using a median smoother on the
#spectral estimate

cost<-function(x,spec,iMax,bLog=TRUE)
{
spec_ar1<-specred(x[1],x[2],spec$freq,0.5/spec$freq[length(spec$freq)])
if (bLog) rmse(log(spec$spec[1:iMax]),log(spec_ar1[1:iMax])) else  rmse(spec$spec[1:iMax],spec_ar1[1:iMax])
}

ar1fit<-function(spec,tss,bLog=TRUE,fmax=NULL)
{
	if (is.null(fmax)) iMax=length(spec$spec) else iMax=which.min(abs(fmax-spec$spec))
	startvalues<-c(sd(tss)^2,acf(tss,plot=FALSE)$acf[1])
	estimated<-nlminb(startvalues,cost,lower=c(0,0),upper=c(1000,0.99),spec=spec,iMax=iMax,bLog=bLog)
	return(estimated)
}



smoothspec<-function(spec,m)
{
	spec$spec=medsmooth(spec$spec,m)$rmedian
	return(spec)
}

specConf<-function(x,p=0.05,p2=0.01,N.R=1000,m=0,spans,fmax=NULL,bLog=F,bMC=F,bMeanspec=F, ...)
{
	spec<-spectrum(x,spans=spans, ...)


	if (m>0)  {
		spec<-smoothspec(spec,m)
		fit<-ar1fit(spec,x,fmax=fmax,bLog=bLog)
	if (bMeanspec) plot(spec,col="red",add=T)
}



	dof<-spec$df
	N.R<-100

if (bMC)
{
	sdsample<-sd(x)
	fixspans=min(length(x)/10,2)
	firstspec<-spectrum(red(fit$par[1],length(x)),spans=fixspans,plot=F)

	specs<-matrix(0,N.R,length(firstspec$spec))
	for (i in 1:N.R) specs[i,]<-spectrum(red(fit$par[1],length(x))*sdsample,plot=F,spans=fixspans)$spec
	meanspec<-apply(specs,2,mean)
} else
{

	meanspec<-specred(fit$par[1],fit$par[2],spec$freq,0.5/spec$freq[length(spec$freq)])

}
	lim<-qchisq(c(1-p/2,p/2),dof)/dof
	lines(spec$freq,meanspec,col="green")
	lines(spec$freq,meanspec*lim[1],col="blue")
	#lines(spec$freq,meanspec*lim[2],col="blue")

	lim<-qchisq(c(1-p2/2,p2/2),dof)/dof
	lines(spec$freq,meanspec,col="green")
	lines(spec$freq,meanspec*lim[1],col="blue",lty=2)
	#lines(spec$freq,meanspec*lim[2],col="blue",lty=2)


}
