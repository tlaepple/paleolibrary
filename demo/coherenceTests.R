### Testcode for coherency.R


result<-vector()
for (i in 2:40) result[i]<-cohbias(0.5,i)
plot(1:40,result,xlab="DOF",ylab="corrected coherence")
abline(h=0.5)

### Testcode for
all<-ts(rnorm(100))
a<-all[-(1:5)]
b<-0*all[-(95:100)]+1*rnorm(95)
result<-    spectrum(cbind(a,b),spans=5)


plot(result,plot.type="coherency")
lines(result$freq,result$coh,col="red",lwd=2)
lines(result$freq,cohbias(sqrt(result$coh)^2,result$df),col="green")

## Confidence limit
abline(h=mcCoherency(a,b,spans=5)^2)

cf<-confPhase(sqrt(result$coh),result$df,debias=FALSE)
cf.debias<-confPhase(sqrt(result$coh),result$df,debias=TRUE)

plot(result,plot.type="phase")
lines(result$freq,result$phase+cf,col="red")
lines(result$freq,result$phase-cf,col="red")
lines(result$freq,result$phase+cf.debias,col="green")
lines(result$freq,result$phase-cf.debias,col="green")

cf.mc<-confPhase.mc(sqrt(result$coh),result$df,spans=5)

plot(cf,type="l",ylim=c(0,5))
lines(cf.mc,col="red")

sum((result$phase+cf)<0)
sum((result$phase+cf.mc)<0)
}


### Now Monte Carlo phase uncertainty

save<-matrix(NA,3,100)
for (i.R in 1:100)
{
 t1<-rnorm(100)
 t2<-rnorm(100)

result<-spectrum(cbind(t1,t2), spans=c(20),log="no",xlim=c(0,0.1))
# plot.spec.phase(result.alf,xlim=c(0,0.1))
save[,i.R]<-confPhase(result,c(1/105,1/41,1/21))$time[,3]
}

hist(save[2,])
hist(save[1,])
hist(save[3,])


### Test sim.coh
## It seems that the bias is not completely removed
spans=21
coh<-c(rep(0,30),rep(1,30))
coh<-abs(sin(1:120/20))
N<-length(coh)*2
s<-rep(0,length(coh))
sb<-rep(0,length(coh))

df<-spectrum(t,spans=spans)$df
for (i.R in 1:100)
{
t<-sim.coh(coh,N)
temp<-spectrum(t,plot.type="coherency",spans=spans)$coh
s<-s+sqrt(temp)
sb<-sb+cohbias(sqrt(temp),df)

}

plot(s/100,type="l",ylim=c(0,1))
lines(sb/100,col="green")
lines(coh,col="red")

test<-sim.coh(1,100)
cor(test[,1],test[,2])

sum<-0
for (i in 1:100)
{
test<-sim.coh(0.5,100)  #Unbiased
sum<-sum+cor(test[,1],test[,2])
}
sum/100



a<-rnorm(100)
b<-abs(fft(a))^2

plot(b[1:50])
lines(b[100:51],col="red")

