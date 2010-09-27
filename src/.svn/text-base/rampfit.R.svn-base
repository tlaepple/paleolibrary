#Implementation of Rampfit
#Ramp (x1,y1 - x2,y2) + two horizontal pieces with a fix length
#Parameters: Search area for x1 and x2, length of horizontal piece

#### Estimates for x1 and x2; Page 3, Mudelsee et al

#sigma2(t) contains the variances
#x(i) the values
#t(i) the time
#i0 = index of the start of flat part of the ramp (in Mudelsee = 1)
#i1  = index of the start of the ramp
#i2  = index of the end of the ramp
#i3   = index of the end of the second flat part (in Mudelsee = n)

ramp_xfit<-function(x,t,sigma2,i0,i1,i2,i3)
{

k1<-sum(1/sigma2[i0:i1])
k2<-sum(1/sigma2[i0:(i2-1)])
k3<-sum(1/sigma2[i2:i3])
k4<-sum(1/sigma2[(i1+1):(i2-1)])
k5<-sum(t[(i1+1):(i2-1)]/sigma2[(i1+1):(i2-1)])
k6<-sum(t[(i1+1):(i2-1)]^2/sigma2[(i1+1):(i2-1)])
k7<-sum(x[i0:(i2-1)]/sigma2[i0:(i2-1)])
k8<-sum(x[i0:i3]/sigma2[i0:i3])
k9<-sum(x[i2:i3]/sigma2[i2:i3])
k10<-sum(x[(i1+1):(i2-1)]/sigma2[(i1+1):(i2-1)])
k11<-sum(x[(i1+1):(i2-1)]*t[(i1+1):(i2-1)]/sigma2[(i1+1):(i2-1)])

t1<-t[i1]
t2<-t[i2]

K1<-k2+(t1*k4-k5)/(t2-t1)  #Equation for the constants from Page 3
K2<-k3-(t1*k4-k5)/(t2-t1)
K3<-k8
K4<-k1+(t2*(t1+t2)*k4+2*k6-(t1+3*t2)*k5)/(t2-t1)^2
K5<-k3+(t1*(t1+t2)*k4+2*k6-(3*t1+t2)*k5)/(t2-t1)^2
K6<-k9-k7-2*(t1*k10-k11)/(t2-t1)

x2<-(K3*K4/K1+K6)/(K2*K4/K1+K5)  #Eq (2)
x1<-(K3-x2*K2)/K1                #Eq (3)

return(list(x1=x1,x2=x2))
}


#Ramp function

#t(i) the time
#i0 = index of the start of flat part of the ramp (in Mudelsee = 1)
#i1  = index of the start of the ramp
#i2  = index of the end of the ramp
#i3   = index of the end of the second flat part (in Mudelsee = n)
#x1  = level of the start of the ramp
#x2  = level of the end of the ramp


ramp<-function(t,i0,i1,i2,i3,x1,x2)
{
	result<-rep(NA,(i3-i0+1))
	result[(i0:(i1-1))-i0+1]<-x1

	result[(i1:i2)-i0+1]<-x1+((t[i1:i2])-t[i1])*(x2-x1)/(t[i2]-t[i1])


	result[((i2+1):i3)-i0+1]<-x2
	return(result)
}




#Rampfit besteht aus einer Brute Force Suche für t und einer LSQ Fit für x (x=Werte)
#x(i) time series
#t(i) the time

#i1_min, i1_max  = index limits of the start of the ramp
#i2_min_i2_max  = index limits of the end of the ramp
#tc1,tc2 = width of the flat part (in time units) 
#sigma2 = uncertainty (variance) for every point in time

rampfit_xt<-function(x,t,i1_min,i1_max,i2_min,i2_max,tc1,tc2,sigma2)
{
	min_mse<-1e9
	for (i1 in i1_min:i1_max)
		for (i2 in i2_min:i2_max)
		{
			i0<-which.min(abs(t-(t[i1]-tc1)))
			i3<-which.min(abs(t-(t[i2]+tc2)))

			if (i0 == 1) warning("Flat part at beginning of the time series")
			if (i3 == length(t)) warning("Flat part at end of the time series")
			
			par<-ramp_xfit(x,t,sigma2,i0,i1,i2,i3)	
		      rampfit<-ramp(t,i0,i1,i2,i3,par$x1,par$x2)
			xpart<-x[i0:i3]
			#plot(t[i0:i3],xpart)
			#lines(t[i0:i3],rampfit)
			mse<-mean((xpart-rampfit)^2)
			if (mse<min_mse)
			{	
				save<-par
				save_i1<-i1
				save_i2<-i2
				min_mse<-mse
				save_i0<-i0
				save_i3<-i3
				
			}
		}
	return(list(i1=save_i1,i2=save_i2,x1=save$x1,x2=save$x2,mse=min_mse,i0=save_i0,i3=save_i3))
}

x<-xClean+rnorm(length(x))*5
result<-rampfit_xt(x,t=seq(x),i1_min,i1_max,i2_min,i2_max,flatwidth=5,sigma2=rep(1,length(x)))


rampfit.val<-function(x,t,param,bPlot=F)
{
	 rampfit<-ramp(t,param$i0,param$i1,param$i2,param$i3,param$x1,param$x2)
	 xpart<-x[param$i0:param$i3]
		
		if (bPlot) {
			plot(t,x);
			lines(t[param$i0:param$i3],rampfit,lwd=2,col="red")
		}
	 return(list(residuals=xpart-rampfit,ramp=rampfit))

}




x<-data
t<-seq(data)
sigma2<-rep(1,length(data))

#Test ramp_xfit(x,t,sigma2,10,31,44,60)  #Geschwindigkeit ca. 4000/Sekunde

plot(ramp(0,10,20,30,5,30))



### 

rampfit<-function(x,t=seq(x),i1_min,i1_max,i2_min,i2_max,tc1,tc2,sigma2=rep(1,length(x)),N.R=10)
{
	bresult<-list()
	bresult$i1<-vector()
	bresult$i2<-vector()
	bresult$x1<-vector()
	bresult$x2<-vector()
	bresult$mse<-vector()


	bestfit<-rampfit_xt(x,t,i1_min,i1_max,i2_min,i2_max,tc1,tc2,sigma2)
	
	bestfit.val<-rampfit.val(x,t,bestfit,bPlot=T)
	a1<-acf(bestfit.val$residuals,plot=F)$acf[2]
	nblock<-lopt(a1,length(bestfit.val$residuals))

	plot(x)
for (i in 1:N.R)
{
	x_sur<-x
	plot(x)

	x_sur[bestfit$i0:bestfit$i3]<-bestfit.val$ramp+blocksample(bestfit.val$residuals,blocklength=nblock) #Surrogates auf alle darauf geben, nicht nur Rampe!
	
	lines(x_sur)

	#bzw... fit auf allen durchführen
	fit.sur<-rampfit_xt(x_sur,t,i1_min,i1_max,i2_min,i2_max,tc1,tc2,sigma2=sigma2)
	bresult$i0[i]<-fit.sur$i0
	bresult$i3[i]<-fit.sur$i3
	bresult$i1[i]<-fit.sur$i1
	bresult$i2[i]<-fit.sur$i2
	bresult$x1[i]<-fit.sur$x1
	bresult$x2[i]<-fit.sur$x2
	bresult$mse[i]<-fit.sur$mse

}
	return(list(bestfit=bestfit,bootstrap=bresult))
}

### 
blocksample<-function(data,blocklength=10)
{
	nblock<-ceiling(length(data)/blocklength)
	starts<-floor(runif(nblock,min=1,max=length(data)-blocklength+1))
	index<-rep(1:blocklength,nblock) + rep(starts,each=blocklength)
	return(data[index[1:length(data)]])
}

lopt<-function(a1,n) {
		if (a1<0) a1=0
		return(max(1,ceiling((6^(0.5)*a1/(1-a1)^2)^(2/3)*n^(1/3))))
		}




## Test code for the ramp
t<-seq(from=1,to=100,by=1)
testramp1<-ramp(t,1,30,60,100,1,10)
sigma2<-rep(1,length(t))

ramp_xfit(testramp1,t,sigma2,1,30,60,100) #OK

#Testcase: non equidistant
index_ne<-cumsum(round(runif(50,min=1,max=5))) #some random sorted indices
testramp_ne<-testramp1[index_ne]  

index<-!is.na(testramp_ne)  #Remove the missing values in x and time
testramp_ne<-testramp_ne[index]
t_ne<-t[index_ne][index]

sigma_ne<-rep(1,length(t_ne)) #No weighting

i0<-1 	#Position of the ramp
i1<-which(testramp_ne>1)[1]-1  
i2<-which(testramp_ne==10)[1]
i3<-length(t_ne)  #

ramp_xfit(testramp_ne,t_ne,sigma_ne,i0,i1,i2,i3) #OK


res<-rampfit_xt(testramp_ne,t_ne,i1_min=1,i1_max=12,i2_min=13,i2_max=19,tc1=10,tc2=10,sigma_ne)


t1<-rampfit.val(testramp_ne,t_ne,res,bPlot=TRUE)

testramp_noise<-testramp_ne+rnorm(length(testramp_ne))
temp<-rampfit(testramp_noise,t_ne,i1_min,i1_max,i2_min,i2_max,tc1,tc2,sigma2=rep(1,length(t)),N.R=100)

plot(t_ne,testramp_noise)
abline(v=t_ne[c(temp$bestfit$i1,temp$bestfit$i2)])

for (i in 1:10)
{
	sr<-ramp(t_ne,temp$bootstrap$i0[i],temp$bootstrap$i1[i],temp$bootstrap$i2[i],temp$bootstrap$i3[i],temp$bootstrap$x1[i],temp$bootstrap$x2[i])
	lines(t_ne[temp$bootstrap$i0[i]:temp$bootstrap$i3[i]],sr,col=i)
}


quantile(temp$bootstrap$x1,c(0.025,0.975))
quantile(temp$bootstrap$x2,c(0.025,0.975))

quantile(temp$bootstrap$i1,c(0.025,0.975))
quantile(temp$bootstrap$i2,c(0.025,0.975))



### Step1: Visual choice of start and endpoints
### Choice of width

## Result: Ramp + uncertainties



index<-!is.na(dO18b)

x<-dO18b[index]*-1
t<-time.huybers[index]

plot(t,x,type="l")
ival<-identify(t,x,n=4)

res<-rampfit(x,t,ival[1],ival[2],ival[3],ival[4],2,2,sigma2=rep(1,length(t)),N.R=10)

abline(v=res$bootstrap$i1,col="red")
abline(v=res$bootstrap$i2,col="blue")

abline(h=res$bootstrap$x1,col="red")
abline(h=res$bootstrap$x2,col="blue")

