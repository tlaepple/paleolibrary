#functions to create surrogate ts


#surrogate timeseries with the same timebasis and autocorrelation as the original
#(based on AR1, / white noise)

snoise.pTs <- function(ts,a1=NULL)
{
if (is.null(a1))
{
ar_result<-ar(ts,order.max=1,aic=FALSE)
nts<-pTs(scale(arima.sim(list(ar = ar_result$ar),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))
}
else nts<-pTs(scale(arima.sim(list(ar = a1),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))

return(nts)
}	

snoise.an.pTs <- function(ts,a1=NULL)
{
if (is.null(a1))
{
ar_result<-ar(ts,aic=TRUE)
nts<-pTs(scale(arima.sim(list(ar = ar_result$ar),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))
}
else nts<-pTs(scale(arima.sim(list(ar = a1),length(ts))),time(ts),getlat(ts),getlon(ts),paste("Surrogate Noise ",getname(ts)))

return(nts)
}	


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
	


#xcorr<- function(x){
#N<-length(x);
#Hxy<-0;
#Rxy<-0;
#for (m in (1:(2*N-1))){
#	if (N-m >=1){
#		Rxy[m]<-0;
#		for (n in (1:(N-m))){
#			Rxy[m]<-(Rxy[m]+(x[n+m]*x[n]))};
#			}
#	if (N-m==0){
#		j<-sum(x*x)}
#	if (m>N){Hxy[m]<-0; Hxy[m]<-Rxy[m-N]}
#	}
#cde<-c(Rxy[length(Rxy):1],j,Hxy[(N+1):(2*N-1)]);
#return(cde/max(cde));
#}
#
#a<-ts(rnorm(5000)+sin(1:5000/5))
#b<-sur.cholesky(a,1)

#spectrum(b,spans=c(3,3))
#spectrum(a,col="red",spans=c(3,3),add=T)

### Normale Korrelation mit kortest; Signifikanz spielt keine Rolle da immer Lokal
#Surrogate 



