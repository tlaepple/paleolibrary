#http://www.rsc.org/images/brief10_tcm18-25920.pdf
#http://finzi.psych.upenn.edu/Rhelp10/2010-February/227865.html
## implemented by Thierry Themenau
# BD Ripley and M Thompson, Regression techniques for the detection
#of analytical bias, Analyst 112:377-383, 1987.

deming<- function(x, y, xstd, ystd, jackknife=TRUE, dfbeta=FALSE,
                   scale=TRUE) {
    Call <- match.call()
    n <- length(x)
    if (length(y) !=n) stop("x and y must be the same length")
    if (length(xstd) != length(ystd)) 
        stop("xstd and ystd must be the same length") 

    # Do missing value processing
    nafun <- get(options()$na.action)
    if (length(xstd)==n) {
        tdata <- nafun(data.frame(x=x, y=y, xstd=xstd, ystd=ystd))
        x <- tdata$x
        y <- tdata$y
        xstd <- tdata$xstd
        ystd <- tdata$ystd
        }
    else {
        tdata <- nafun(data.frame(x=x, y=y))
        x <- tdata$x
        y <- tdata$y
        if (length(xstd) !=2) stop("Wrong length for std specification")
        xstd <- xstd[1] + xstd[2]*x
        ystd <- ystd[1] + ystd[2] * y
        }

    if (any(xstd <=0) || any(ystd <=0)) stop("Std must be positive")

    minfun <- function(beta, x, y, xv, yv) {
        w <- 1/(yv + beta^2*xv)
        alphahat <- sum(w * (y - beta*x))/ sum(w)
        sum(w*(y-(alphahat + beta*x))^2)
        }

    minfun0 <- function(beta, x, y, xv, yv) {
        w <- 1/(yv + beta^2*xv)
        alphahat <- 0  #constrain to zero
        sum(w*(y-(alphahat + beta*x))^2)
        }

    afun <-function(beta, x, y, xv, yv) {
        w <- 1/(yv + beta^2*xv)
        sum(w * (y - beta*x))/ sum(w)
        }

    fit <- optimize(minfun, c(.1, 10), x=x, y=y, xv=xstd^2, yv=ystd^2)
    coef = c(intercept=afun(fit$minimum, x, y, xstd^2, ystd^2), 
               slope=fit$minimum)
    fit0 <- optimize(minfun0, coef[2]*c(.5, 1.5), x=x, y=y, 
                     xv=xstd^2, yv=ystd^2)

    w <- 1/(ystd^2 + (coef[2]*xstd)^2) #weights
    u <- w*(ystd^2*x + xstd^2*coef[2]*(y-coef[1])) #imputed "true" value
    if (is.logical(scale) && scale) {
        err1 <- (x-u)/ xstd
        err2 <- (y - (coef[1] + coef[2]*u))/ystd
        sigma <- sum(err1^2 + err2^2)/(n-2)
        # Ripley's paper has err = [y - (a + b*x)] * sqrt(w); gives the same SS
        }
    else sigma <- scale^2
    
    test1 <- (coef[2] -1)*sqrt(sum(w *(x-u)^2)/sigma) #test for beta=1
    test2 <- coef[1]*sqrt(sum(w*x^2)/sum(w*(x-u)^2) /sigma) #test for a=0
                      
    rlist <- list(coefficient=coef, test1=test1, test0=test2, scale=sigma,
                  err1=err1, err2=err2, u=u)

    if (jackknife) {
        delta <- matrix(0., nrow=n, ncol=2)
        for (i in 1:n) {
            fit <- optimize(minfun, c(.5, 1.5)*coef[2], 
                            x=x[-i], y=y[-i], xv=xstd[-i]^2, yv=ystd[-i]^2)
            ahat <- afun(fit$minimum, x[-i], y[-i], xstd[-i]^2, ystd[-i]^2)
            delta[i,] <- coef - c(ahat, fit$minimum)
            }
        rlist$variance <- t(delta) %*% delta
        if (dfbeta) rlist$dfbeta <- delta
        }

    rlist$call <- Call
    class(rlist) <- 'deming'
    rlist
    }


### Monte Carlo Comparison of WLSQ and OLSBC
#r1<-r2<-vector()
#for (i in 1:100)
#{
#sx<-50
#sy<-10
#x<-1:100+rnorm(100)*sx
#y<-1:100+rnorm(100)*sy

#kappa<-1/((1-sx^2/var(x)))
#r1[i]<-lm(y~x)$coeff[2]*kappa
#r2[i]<-deming(x,y,rep(sx,100),rep(sy,100))$coef[2]
#}

