#Library for multiple plotting and picking

#Input; List of X and List of Y, xlim can be a list or a single xlim range
PlotMultiple<-function(x,y,type="b",xlab="time",xlim=range(x),bPoints=TRUE,pch=20,...)
{
	windows()
	close.screen(1,all.screens = TRUE)
	if (!is.list(y)) stop("y has to be a list")
                xlimList<-NULL
                if (is.list(xlim)) xlimList<-xlim  #multiple xlims
	N<-length(y)
	if (!is.list(x))
		{
			xnew<-list()
			for (i in 1:N) xnew[[i]]<-x
			x<-xnew
		}

	split.screen(c(N+1,1))

	for (i in 1:N)
	{
		screen(i)
		par(mai=c(0,0,0,0.1))
                                if (!is.null(xlimList)) xlim<-xlimList[[i]]
		if (i<N) {
                      plot(x[[i]],y[[i]],axes=F,type=type,xlim=xlim,...)
                      if (bPoints) {
                          points(x[[i]],y[[i]],col="grey",pch=pch)
                          lines(x[[i]],y[[i]],...)
                      }
                       } else
	{
		plot(x[[i]],y[[i]],axes=F,type=type,xlim=xlim,xlab="time kyr BP",...)
                if (bPoints) {
                    points(x[[i]],y[[i]],col="grey",pch=pch)
                     lines(x[[i]],y[[i]],...)
                }
			axis(1)
			}
	box()
	}


}

#Identifies Points and gives them colors
identifyPch <- function(x, y=NULL, n=6, pch=19,colors=rep(c("blue","green","red"),2), ...)
{

    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    NId=0;
 res <- integer(0)

    while(NId < n) {
        ans <- identify(x, y, n=1, plot=FALSE, ...)
        if(!length(ans)) break
        points(x[ans], y[ans], pch = pch,col=colors[NId+1])
        res <- c(res, ans)
	  NId<-NId+1

    }
    res
}



# Pick NPick values at row i
PickN<-function(x,y,i,NPick=6)
{
	if (!is.list(x))
		{
			xnew<-list()
			for (i in 1:N) xnew[[i]]<-x
			x<-xnew
		}



	screen(i)
	return(identifyPch(x[[i]],y[[i]],n=NPick))
}

