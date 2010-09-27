#Implementation of standard finite response filters, e.g. Bloomfield 1976, linear filtering

#Derives and plots the transfer function (given a filter)
get.transfer<-function(g.u,resolution=100,plot=T,add=F, ...)
{
#Get the transfer function of a symetric filter, page 122
#n<-length of the filter (odd number)
n<-length(g.u)
n.side<-(n-1)/2
omega=(1:resolution)*pi/resolution
yt<-rep(0,length(omega))

for (u in (-1*n.side):n.side)
	yt<-yt+(g.u[u+n.side+1]*exp(-1*1i*omega*u)) 

if (plot) {
	if (!add) plot(omega/2/pi,abs(yt)^2,type="l",xlab="omega",main="Transfer function", ...)
		else lines(omega/2/pi,abs(yt)^2,col="blue", ...)
}

return(list(omega=omega/2/pi,y=abs(yt)^2))
}

#Derive the (smoothed) least square lowpass, given the cutoff frequency omega.c and the length of the filter n
lowpass<-function(omega.c,n=9,sample=1,convergence=T)
{
if ((n %% 2) == 0) stop("N must be odd, this function calculates only symetrical = phase preserving filters")

	omega.c<-omega.c/sample
if (omega.c >= 0.5) stop("frequency higher or equal then Nyquist")
#calculate least square solution for ideal lowpass
	g.u<-0
	n.side<-(n-1)/2
	for (u in 1:n.side)
	{
		g.u[u+n.side+1]<-sin(u*omega.c*2*pi)/(pi*u)
	}
	#g_0
	g.u[n.side+1]<-omega.c*2
	#Now flip over
	g.u[1:n.side]<-g.u[n:(n.side+2)]

if (convergence)
{
	#Multiply by convergence factors to reduce the ripples
	delta<-4*pi/n
	for (u in c((-1*n.side):-1,1:n.side)) g.u[u+n.side+1]<- g.u[u+n.side+1]*sin(u*delta/2)/(u*delta/2)
}	
	
	return(g.u)
}


#Derive the smoothed least square highpass, given the cutoff frequency omega
highpass<-function(omega.c,n=9,sample=1,convergence=T)
{
	directtransfer<-rep(0,n)
	directtransfer[(n+1)/2]<-1
	return(directtransfer-lowpass(omega.c,n=n,sample=sample,convergence=convergence))
}




bandpass<-function(omega.upper,omega.lower,n=9,sample=1,convergence=T)
{
if ((n %% 2) == 0) stop("N must be odd, this function calculates only symetrical = phase preserving filters")

	
	return(lowpass(omega.upper,n,sample,convergence)-lowpass(omega.lower,n,sample,convergence))
}

filter.pTs <- function(data,filter,...)
{
	result<-filter(data,filter,...)
	return(pTs(result,time(data),getlat(data),getlon(data),paste(getname(data),"filtered"),gethistory(data)))
}

#Endpoint problem e.g. Mann et al., GRL 2003
#While constrain ts (1) – (3) can be applied explicitly in
#the frequency domain [ e.g., Park , 1992; Ghil et al. , 2002], it
#is possible to imple ment reasonable approximations to these
#co nstraints in a simple manner in the time domain as
#follows: To approximate the ‘minimum norm’ constraint,
#one pads the series with the long-term mean beyond the
#boundaries (up to at least one filter width) prior to smoothing. 
#To approximate the ‘minimum slope’ constraint, onepads the series 
#with the values within one filter width of theboundary reflected about the time boundary.
#This leads thesmooth towards zero slope as it approaches the boundary.
#Finally, to appr oxima te the ‘ minimum r oughness’ constraint,
# one pads the series with the values within one filterwidth of the boundary 
#reflected abou t the time boundary,and reflected vertically (i.e., about the ‘‘y’’ axis) 
#relativ e tothe final value. This tends to impose a point of inflection atthe boundary, 
#and leads the smooth towards the boundarywith constant slope.

#Here write the filter with the modified boundary conditions
filter.pTs1 <- function(data,filter,method=1,...)
{
N<-floor(length(filter)/2)
	if (method == 1)  #Minimum Norm
	{
		before<-rep(mean(data),N)
		after<-rep(mean(data),N)
	}
	if (method == 2)
	{
		before<-c(data)[N:1]
		after<- c(data)[length(data):(length(data)-N+1)]
	}
	if (method == 3)
	{
		before<-c(data)[N:1]
		after<- c(data)[length(data):(length(data)-N+1)]

		before<-c(data)[1]-(before - mean(before))
		after<-c(data)[length(data)]-(after - mean(after))



	}
	
	result<-filter(c(before,data,after),filter,circular=F)[(N+1):(N+length(data))]
	return(pTs(result,time(data),getlat(data),getlon(data),paste(getname(data),"filtered"),gethistory(data)))
}


SSA2 <- function(ts,L,I=1,plot=F) 
{
x1<-na.omit(as.vector(ts))
#create the trajectory matrix
N=length(x1)
 if (L>N/2) L<-N-L
K<-N-L+1
X = matrix(0,L,K)
for (i in 1:K) X[1:L,i]<-x1[i:(L+i-1)]

S<-X%*%t(X)
 ev<-eigen(S)
U<-ev$vectors
d<-diag(ev$values)
sev<-sum(ev$values)
 if (plot) plot((ev$values/sev)*100)
V<-t(X)%*%U
rc<-U%*%t(V)

Vt<-as.matrix(t(V)) 
rca<-as.matrix(U[,I])%*%Vt[I,] 
y<-rep(0,N)
Lp<-min(L,K) 
Kp<-max(L,K)
for (k in (0:(Lp-2)) )
{
	 for (m in (1:(k+1)) ) y[k+1]<-y[k+1]+(1/(k+1))*rca[m,(k-m+2)]
} 
for (k in ((Lp-1):(Kp-1)) )
{
	for (m in (1:Lp))  y[k+1]<-y[k+1]+(1/(Lp))*rca[m,(k-m+2)];
}
dm<-dim(rca)
for ( k in (Kp:N))
{
	 for (m in ((k-Kp+2):(N-Kp+1)))  {
		if ((m <=dm[1]) && ((k-m+2)<=dm[2]))  y[k+1]<-y[k+1]+(1/(N-k))*rca[m,(k-m+2)]
} }

ts_result<-pTs(y,time(ts),0,0,getname(ts))
invisible(ts_result)

 }



