
# Plot ECDF of example from notes
ee <- ecdf(c(5,3,7))
plot(ee)

# Now plot ECDF for nerve data from text
x <- scan("http://www.stat.cmu.edu/~larry/all-of-statistics/=data/nerve.dat")
ee <- ecdf(x)
plot(ee)

ee(.5)
abline(v=.5)
abline(h=ee(.5))

ee(.6) - ee(.4)

# Draw our own and add nonparametric confidence band. Adapted from R
# code on author's website:
# http://www.stat.cmu.edu/~larry/all-of-statistics/=Rprograms/edf.r
ecdf2 = function(x,CI=TRUE) {
  ox <- sort(x) # ordered x
  n <- length(x)
  Fhatn <- (1:n)/n
  plot(ox,Fhatn,type="l")
  rug(x,ticksize=0.025)
  if(CI) {
    alpha <- 0.05
    eps <- sqrt(log(2/alpha)/(2*n))
    upper <- pmin(Fhatn + eps,1)
    lower <- pmax(Fhatn - eps,0)
    lines(ox,upper,type="s",lwd=2,col=2,lty=2)
    lines(ox,lower,type="s",lwd=2,col=2,lty=2)
  }
}
ecdf2(x)

# Simulate some uniform(0,1) data and draw sample paths.
set.seed(8675309)
n=10
ecdf2(runif(n))
abline(a=0,b=1,lwd=2,col="blue") # add true CDF
n=100; ecdf2(runif(n)); abline(a=0,b=1,lwd=2,col="blue") 
n=1000; ecdf2(runif(n)); abline(a=0,b=1,lwd=2,col="blue") 
n=10000; ecdf2(runif(n)); abline(a=0,b=1,lwd=2,col="blue") 
# ECDF getting closer to true CDF as n increases. 

# Centre ECDFs by true CDF to focus on differences.
c.ecdf = function(x,F) {
  ox <- sort(x)
  n <- length(x)
  Fhatn <- (1:n)/n
  ce = Fhatn - F(ox)
  plot(ox,ce,type="l",ylim=c(-.2,.2))
  rug(x,ticksize=0.025)
}

n=10; c.ecdf(runif(n),punif)
n=100; c.ecdf(runif(n),punif)
n=1000; c.ecdf(runif(n),punif)
n=10000; c.ecdf(runif(n),punif)
n=100000; c.ecdf(runif(n),punif)
# ECDF converging to CDF -- Glivenko Cantelli says this convergence
# is uniform and almost sure.

# Scale up Fhat-F by sqrt(n); sqrt(n)(Fhat-F) is called the 
# empirical process.
c.ecdf = function(x,F) {
  ox <- sort(x)
  n <- length(x)
  Fhatn <- (1:n)/n
  ce = sqrt(n) *(Fhatn - F(ox))
  plot(ox,ce,type="l",ylim=c(-1,1))
  rug(x,ticksize=0.025)
}

n=10; c.ecdf(runif(n),punif)
n=100; c.ecdf(runif(n),punif)
n=1000; c.ecdf(runif(n),punif)
n=10000; c.ecdf(runif(n),punif)
n=100000; c.ecdf(runif(n),punif)

# Donsker's Theorem says that the empirical process "converges in distribution"
# to a Gaussian process called a Brownian bridge, or "tied down"
# Brownian motion. (Notice how the process is always 0 at x=0 and 1.)
# What does convergence of stochastic processes mean? 