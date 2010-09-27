#scale each timeseries=point of the field

#Here only the history gets modified by removing the last 1 / 2 entries and
#added by Ops.pField and adding a new one... 
scale.pField <- function(x, center = TRUE, scale = TRUE)
{
    result<-NextMethod()
    hist<-attr(result,"history")
    hist<-hist[1:(length(hist)-sum(center,scale))]
    hist<-c(hist,paste(date(),"scale center=",center,"scale=",scale))
    attr(result,"history")<-hist
    return(result)
}


scale.pTs <- function(x, center = TRUE, scale = TRUE)
{
    hist<-gethistory(x)
    hist<-hist[1:(length(hist)-1)]
    result<-pTs(NextMethod(),time(x),getlat(x),getlon(x),getname(x),hist,date=FALSE)
    hist<-c(hist,paste(date(),"scale center=",center,"scale=",scale)) 
    return(addhistory(result,hist))
}



detrend <- function(x, ...) UseMethod("detrend")

detrend.default <- function(x, ...) print("No default detrend function")



detrend.pField<- function(x)
{
	x[]<-lm(unclass(x)~seq(nrow(x)))$residuals
	return(addhistory(x,"detrend"))
}



detrend.pTs<- function(x)
{
	x[]<-lm(unclass(x)~seq(length(x)))$residuals
	return(addhistory(x,"detrend"))

}



rollmean.pField <- function(x,k,na.pad = FALSE, ...)
{
  x[]<-rollmean.default(x,k,na.pad=TRUE, ...)
  if (na.pad==FALSE) return(addhistory(na.omit(x),paste("rollmean",k))) 
	else return(addhistory(x,paste("rollmean",k)))

}


rollmean.pTs <- function(x,k,na.pad = FALSE, ...)
{
  x[]<-rollmean.default(x,k,na.pad=TRUE, ...)
  if (na.pad==FALSE) return(addhistory(na.omit(x),paste("rollmean",k))) 
	else return(addhistory(x,paste("rollmean",k)))

}



fastcor.pTs <- function(pTs,pField)
{
return(pField(cor(pTs,pField),9999,getlat(pField),getlon(pField),paste("cor",getname(pTs),getname(pField)),gethistory(pField),date=FALSE))
}

#Correlaties pTs with pField or pTs with pTs
cor.pTs <- function(pTs,pField,debug=F)
{
	#bring both on the same time basis
	start<-max(start(pTs)[1],start(pField)[1])
	end<-min(end(pTs)[1],end(pField)[1])
	if (debug) print(paste("Common time period: ",start,end))
	pTs<-window(pTs,start,end)
	pField<-window(pField,start,end)

	if (class(pField)[1]=="pField") #field correlation
	{

		class(pField)<-"matrix"

		#Filter out data with contain only missing values
		dat<-pField[,!is.na(colSums(pField))]
		result<-matrix(NA,1,ncol(pField))
	
		tresult<-cor(dat,pTs)
		result[,!is.na(colSums(pField))]<-tresult
		class(pField)<-"pField"
		return(copyattr(result,pField))
	}
	else return(cor(pTs,pField))
}


cortest.pTs <- function(pTs,pField)
{
#bring both on the same time basis
start<-max(start(pTs)[1],start(pField)[1])
end<-min(end(pTs)[1],end(pField)[1])
print(paste("Common time period: ",start,end))
pTs<-window(pTs,start,end)
pField<-window(pField,start,end)
class(pField)<-"matrix"

#Filter out data with contain only missing values
dat<-pField[,!is.na(colSums(pField))]
result<-matrix(NA,2,ncol(pField))

tresult<-apply(dat,2,mycor.test,c(pTs))
result[,!is.na(colSums(pField))]<-tresult


class(pField)<-"pField"
return(copyattr(result,pField))
}



mycor.test <- function(v1,v2, ...) 
{
	t<-cor.test(v1,v2, ...) 
	return(c(t$estimate,t$p.value))
}

prcompO.pTs <- function(data,nPc=2,center=TRUE,scale=TRUE, ...)
  {
      temp<-attributes(data)
        class(data)<-"matrix"
        result<-prcomp(data,center=center,scale=scale)
        pc<-pTs(result$x[,1:nPc],time(data),9999,9999,paste("PC",1:nPc),gethistory(data),date=FALSE)
        pc<-addhistory(pc,"PRCOMP")
        eof<-result$rotation[,1:nPc]
        var<-result$sdev[1:nPc]^2
        varsum<-sum(result$sdev^2)
        return(list(pc=pc,eof=eof,var=var,varsum=varsum))
    }


prcompO.pField <- function(data,nPc=10,center=TRUE,scale=TRUE, ...)
{
  temp<-attributes(data)
  class(data)<-"matrix"
  result<-prcomp(data,center=center,scale=scale)
  pc<-pTs(result$x[,1:nPc],time(data),9999,9999,paste("PC",1:nPc,temp$name),gethistory(data),date=FALSE)
  pc<-addhistory(pc,"PRCOMP")
  eof<-pField(result$rotation[,1:nPc],1:nPc,temp$lat,temp$lon,paste("EOF",temp$name),gethistory(data),date=FALSE)
  eof<-addhistory(eof,"PRCOMP")
  var<-result$sdev[1:nPc]^2
  varsum<-sum(result$sdev^2)
  return(list(pc=pc,eof=eof,var=var,varsum=varsum))
}




composite.pTs <- function(ts,field,sign=FALSE,pval=TRUE,sp=sd(ts),sm=(-1*sp))
{
	temp<-attributes(field)
      field<-scale(field,scale=FALSE)
	ts<-scale(ts)
#bring both on the same time basis
	start<-max(start(ts)[1],start(field)[1])
	end<-min(end(ts)[1],end(field)[1])
	print(paste("Common time period: ",start,end))
	ts<-window(ts,start,end)
	field<-unclass(window(field,start,end))

	ts[is.na(ts)]<-0
	if (sign)
	{
		index.plus<-ts>0
		index.minus<-ts<0
	}
	else
	{
		
		
		index.plus <- ts>sp
		index.minus <- ts<sm
	}
	
	field.plus<-field[index.plus,]
	field.minus<-field[index.minus,]
	if (pval)
	{
	field.p.plus<-rep(NA,ncol(field))
	for (i in 1:ncol(field)) field.p.plus[i]<-t.test(field.plus[,i])$p.value
	

        field.p.minus<-rep(NA,ncol(field))
	for (i in 1:ncol(field)) field.p.minus[i]<-t.test(field.minus[,i])$p.value

 	
      field.p.symetrie<-rep(NA,ncol(field))
	for (i in 1:ncol(field)) field.p.symetrie[i]<-t.test(-1*field.minus[,i],field.plus[,i])$p.value	
         }
        

        field.plus<-pField(colMeans(field.plus),9999,temp$lat,temp$lon,paste("+ Compos.",getname(ts),temp$name),temp$history,date=FALSE)
	field.minus<-pField(colMeans(field.minus),9999,temp$lat,temp$lon,paste("- Compos.",getname(ts),temp$name),temp$history,date=FALSE)

        field.p.plus<-pField(field.p.plus,9999,temp$lat,temp$lon,paste("P Compos.+",getname(ts),temp$name),temp$history,date=FALSE)
	field.p.minus<-pField(field.p.minus,9999,temp$lat,temp$lon,paste("P Compos.-",getname(ts),temp$name),temp$history,date=FALSE)
	field.p.symetrie<-pField(field.p.symetrie,9999,temp$lat,temp$lon,paste("P Compos.-",getname(ts),temp$name),temp$history,date=FALSE)
	
	return(list(plus=field.plus,minus=field.minus,p.plus=field.p.plus,p.minus=field.p.minus,p.symetrie=field.p.symetrie))
	

}


pcor <- function(data,p1,p2,t1,t2)
{
	data<-scale(detrend(window(data,t1,t2)))
	x<-cortest.pTs(selspace(data,lat1=p1,lon1=p2),data)
	result<-x[1,]
	result[x[2,]>0.05]<-NA
	return(result)
}

sigcor <- function(ts,data,p=0.05)
{
	
	x<-cortest.pTs(ts,data)
	result<-x[1,]
	result[x[2,]>p]<-NA
	return(result)
}

####Experimental unfinished functions

ccf.pTs <- function(ts,field)
{
	result<-pField(NULL,9999,getlat(field),getlon(field),"CCF ANTICOR LAG")
	for (i in 1:ncol(field)) 
		result[i]<-myccf(field[,i],ts)
	return(result)
}
myccf<-function(ts1,ts2)
{
   temp<-ccf(ts1,ts2,plot=FALSE,lag.max=14)
   return(temp$lag[temp$acf==min(temp$acf)][1])
}


ccf.pTs<-function(pTs1,pTs2,debug=F,...)
{
	#bring both on the same time basis
	start<-max(start(pTs1)[1],start(pTs2)[1])
	end<-min(end(pTs1)[1],end(pTs2)[1])
	if (debug) print(paste("Common time period: ",start,end))
	pTs1<-window(pTs1,start,end)
	pTs2<-window(pTs2,start,end)
	return(ccf(c(pTs1),c(pTs2)))
}

### Some Pearson correlation stuff
#Both sided p-value for Pearson correlation
pval.cor<-function(r,N)
{
df<-N-2
t=r*sqrt(df)/sqrt(1-r^2)
return((1-pt(t, df, lower.tail = T) )+pt(t, df, lower.tail = F))
}

#Given the p-value (both sided), and N
# in returns r for Pearson correlation
rval.cor<-function(p,N)
{
df<-N-2
t<-qt(p/2,df,lower.tail=F)
r<-(0:999)/1000
return(r[which.min(abs(t-r*sqrt(df)/sqrt(1-r^2)))])
}



