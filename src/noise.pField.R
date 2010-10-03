#functions to create surrogate ts

#surrogate timeseries with the same timebasis and autocorrelation (AR1) as the original
#if order.max > 1, or order.max =NULL, than an AR-N process is fitted
#if a1 is given, an AR1 process with a1 is created.
#in the moment, the result is scaled (sd=1,mean=0)
#ToDo: Add variance scaling

snoise.pTs <- function(ts,a1=NULL,order.max=1)
{
if (is.null(a1))
{
if (order.max != 1) ar_result<-ar(ts,order.max=order.max,aic=TRUE) else ar_result<-ar(ts,order.max=1,aic=FALSE)
nts<-pTs(scale(arima.sim(list(ar = ar_result$ar),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))
}
else nts<-pTs(scale(arima.sim(list(ar = a1),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))

return(nts)
}



### return ACF, timelag 1
get.a1<-function(x) return(acf(c(x))$acf[2])

### simulate a realization of an AR1 process with a1
red<-function(a1,n) return(c(arima.sim(list(ar = a1),n)))

#create noise with the same distribution, by folding a uniform distribution with
#the quantile function
#up to now only for independent samples
#only for univariate ts
enoise.pTs <- function(ts)
{
n<-length(time(ts))
noise_data<-ts
q<-quantile(c(noise_data),probs=(1:1000)/1000)
noise_sim<-q[runif(n)*1000+1]
	#correct the noise to get the same sd
noise_sim<-noise_sim*(sd(noise_data)/sd(noise_sim))
return(pTs(noise_sim,time(ts),name="enoise"))
}



#create N.R random timeseries with the same autocovariance matrix as ts_in
#see Haam and Huybers, PO 2010
sur.cholesky<-function(ts_in, N.R){

### function to generate N.R realizations.

		ts_work<-scale(ts_in)
		TS<-matrix(0, nrow=length(ts_work), ncol=N.R);
		cv<-acf(ts_work,lag.max=length(ts_work),plot=F)$acf;

	#Build autocovariance matrix !?
		cvcs<-matrix(0, nrow=length(cv), ncol=length(cv));
		for (i in 1:(length(cv)-1)) cvcs[,i]<-c(cv[(i:1)],cv[2:(length(cv)-i+1)])
	 	cvcs[,length(cv)]<-rev(cv)



			#erster Part: Umgedreht: 1tes, 2.1. 3.2.1. 4.3.2.1.
			# dann alle bis auf die letzten i



		U<-chol(cvcs);	# cholesky decomposition.

		M<-matrix(rnorm(N.R*max(dim(U))), nrow=N.R, ncol=max(dim(U)));
		TS<-t(M%*%U);


}


###Simulate two timeseries, with length N, and the coherency cb
## Coherency has to be a scalar, or a vector with half the length of N containing the
##Coherency against frequency
## Description in Huybers, Nature Geoscience 2008, Supplement, or in Wunsch timeseries primer
sim.coh<-function(cb,N=1000)
{
if (length(cb)==1) cb<-rep(cb,N/2)
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



### Testcode for the AR1 processes

#test<-red(0.9,1000)
#get.a1(test)
#get.a1(snoise.pTs(test))
#get.a1(red(get.a1(test),1000))

# Tescode for sur.cholesky
#a<-ts(rnorm(5000)+sin(1:5000/5))
#b<-sur.cholesky(a,1)

#spectrum(b,spans=c(3,3))
#spectrum(a,col="red",spans=c(3,3),add=T)

### Normale Korrelation mit kortest; Signifikanz spielt keine Rolle da immer Lokal
#Surrogate



