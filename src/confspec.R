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

#x2 gets replaced by noise
mcCoherency<-function(x1,x2,spans=2,N.R=100,p=0.95)
{
    test<-spectrum(x1)
    cohsave<-matrix(NA,N.R,length(test$freq))
    for (i.R in 1:N.R)
    {
        temp<-spectrum(cbind(snoise.pTs(x2),x1),spans=spans,plot=F)
       # plot.spec.coherency(temp)
        cohsave[i.R,]<-sqrt(temp$coh)
      }
   return(quantile(cohsave,p))
}

mcCoherency(all.lr04[,5],all.lr04[,1],spans=c(10))
mcCoherency(w(alf.lr04),w(all.lr04[,1]),spans=c(10))

iNext<-function(cx,x)
{
    return(which.min(abs(cx-x)))
}

confPhase<-function(x,freq,ci=0.95)
{
   gg <- 2/x$df
        coh <- sqrt(x$coh)
        cl <- asin(pmin(0.9999, qt(ci, 2/gg - 2) * sqrt(gg *
            (coh^{
                -2
            } - 1)/(2 * (1 - gg)))))


result.angle<-matrix(NA,3,length(freq))
result.time<-result.angle
   result.coh<-vector()
   for (i in 1:length(freq))
   {
       index<-iNext(x$freq,freq[i])
       result.angle[1,i]<-x$phase[index]-cl[index]
       result.angle[2,i]<-x$phase[index]
       result.angle[3,i]<-x$phase[index]+cl[index]
       result.time[,i]<-result.angle[,i]/2/pi*(1/freq[i])
       result.coh[i]<-sqrt(x$coh[index])
   }
   return(list(angle=result.angle,time=result.time,coh=result.coh))
}

### Validieren
#all<-ts(rnorm(1000))
#a<-all[-(1:10)]
#b<-all[-(991:1000)]
 #  result<-    spectrum(cbind(a,b),spans=10)
#plot(result,plot.type="phase")
#confPhase(result,0.04)

save<-matrix(NA,3,100)
for (i.R in 1:100)
{
 t1<-pTs(rnorm(length(w(alf.lr04))),time(w(alf.lr04)))
 t2<-0.5*t1+pTs(rnorm(length(w(alf.lr04))),time(w(alf.lr04)))

result.alf<-spectrum(cbind(t1,t2), spans=c(20),log="no",xlim=c(0,0.1))
# plot.spec.phase(result.alf,xlim=c(0,0.1))
save[,i.R]<-confPhase(result.alf,c(1/105,1/41,1/21))$time[,3]
}

hist(save[2,])
hist(save[1,])
hist(save[3,])

sum(save[3,]<0)

plot(result.alf$coh)


###Simulate two timeseries, with length N, and the coherency cb
##
simWCoherency<-function(cb,N=1000)
{

    if (length(cb)!=(N/2)) stop("coherency has to have the length N/2")
    cb[(N+1):(2*N)]<-cb[N:1]


  fx<-fft(rnorm(N));
  fx<-fx/sum(abs(fx));
  fy<-fft(rnorm(N));
  fy<-fy/sum(abs(fy));
  ys <-Re(fft(fy*sqrt(1-Conj(cb)^2),inverse=TRUE))/length(fy)
  ys =ys+Re((fft((fx*Conj(cb)),inverse=TRUE)))/length(fx);
  xs =Re((fft(fx,inverse=TRUE))/length(fx));
 return(cbind(xs,ys))
}


### Read the coherency bias mat file

library(R.matlab)
cohraw<-readMat("e:/data/paleoLibrary/data/cohbias.mat")
save(cohraw,file="e:/data/paleoLibrary/data/cohbias.dat")

cb<-0.5
v<-0.8

### Interpoliere die ensprechenden Stellen in y Richtung zum richtigen Freiheitsgrad

ec<-vector()
### Interpolieren zum richtigen Freiheitsgrad
for (i in 1:length(cohraw$c)) ec[i]<-approx(cohraw$n,cohraw$expect[,i],v)$y
#Interpolieren zu cb
approx(ec,cohraw$c,cb)$y

end;


cu=cu(:);

pl=find(isnan(cu)==1 & cb<1 & cb>=0); %If cu is NaN while cb is between (0,1)
cu(pl)=0;


spectrum(rnorm(1000),spans=10)$df
