#Small shortcuts
last<-function (x)
{
    return(x[length(x)])
}

first<-function (x)
{
    return(x[1])
}

#Search the element in xvector which is closest to x
#if N = nearest; L = less, the next <= smaller value; M=more, the >= next higher value
closest.element<-function(xvector,x,type="N")
{
   if (type=="N")  return(which.min(abs(x-xvector)))

   if (min(diff(xvector))<0) stop("Vector must be monotonically increasing for the the methods M and L")
   if (type=="M") return(first(which(x<=xvector)))
   if (type=="L") return(last(which(x>=xvector)))

}

binAvg<-function(x,y,N=NULL,breaks=pretty(x,N),bFill=FALSE)

#Averages y into bins according to the positon of a in the breaks
#Either give N=number of breaks, or N+1 breaks
#Breaks are defined as x>breaks[i], and x<=breaks[i+1]
#see also stats.bin in library(fields) which returns a more comprehensive statistics
#if fill=T, fill empty bins using linear interpolation from the neighbours to the center of the bin
#Returns the breaks, centers, the averaged values and nobs, the number of observations averages
{
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    avg<-rep(NA,length(breaks)-1)
    nobs<-rep(NA,length(breaks)-1)
    for (i in 1:(length(breaks)-1)) {
        selection<-y[which((x>breaks[i])&(x<=breaks[i+1]))]
        avg[i]<-mean(na.omit(selection))
        nobs[i]<-sum(!is.na(selection))
        }


  if ((sum(is.na(avg))>0)&(bFill))
  {
      yInt<-approx(x,y,centers)$y
      missing<-is.na(avg)
      avg[missing]<-yInt[missing]
  }
   return(list(breaks=breaks,centers=centers,avg=avg,nobs=nobs))
}

#SINC function


sinc<-function (X)
{
    Y <- X
    Z <- X == 0
    Y[Z] <- 1
    Y[!Z] <- sin(pi * X[!Z])/(pi * X[!Z])
    return(Y)
}


### Fill missing values using linear interpolation
na.fill<-function(x) {
    if (sum(is.na(x))==0) return(x)
    nonMissing<-!is.na(x)
     return(pTs(approx(c(time(x))[nonMissing],c(x)[nonMissing],c(time(x)),rule=2)$y,time(x),getlat(x),getlon(x)))
               }


##Source a whole directory
 sourceDir <- function(path, trace = TRUE, ...) {
         for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
            if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
         }
      }

### approximate a timeseries using the nearest neighbour
## just extends approx which always takes the right or left neighbour .. or the weighted mean between both if f>0<1
approx.nearest<-function(x,y,xi)
{
result<-list()
result$x=xi
result$y= approx(c(x[1], x + c(diff(x)/2, 0)), c(y[1], y), xi,
            method = "constant", f = 1)$y
return(result)
    }


unfactor<-function(f)
#Gets the numeric value from a factor
    #see ?factor
    {
       return(as.numeric(levels(f)))
   }
