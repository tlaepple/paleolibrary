

binSpec<-binAvg(1:10,1:10,breaks=0:100,bFill=T)
nBins<-100
spec<-spectrum(cdat[[1]][,100])
plot(spec$freq,spec$spec,type="l",log="xy")
 binSpec<-binAvg(log(spec$freq),spec$spec,breaks=breaks,bFill=T)
 index<-!is.na(binSpec$avg)
    binSpec$avg=binSpec$avg[index]
    binSpec$centers=binSpec$centers[index]
    binFreq<-exp(binSpec$centers)  #Get back of the

    plot(spec$freq,spec$spec,log="xy",type="l",xlab="f (1/kyr)",ylab="PSD")
    lines(binFreq,binSpec$avg,col="red")
abline(v=c(0.8,1.2),col="green")
#abline(v=exp(binSpec$center))


breaks=seq(from=-7.85,to=1.75,by=0.1)
