
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
ox <- sort(x) # ordered x
n <- length(x)
Fhatn <- (1:n)/n
plot(odat,Fn,type="l")
alpha <- 0.05
eps <- sqrt(log(2/alpha)/(2*n))
print(eps)
upper <- pmin(Fhatn + eps,1)
lower <- pmax(Fhatn - eps,0)
lines(ox,upper,type="s",lwd=2,col=2,lty=2)
lines(ox,lower,type="s",lwd=2,col=2,lty=2)
rug(x,ticksize=.025)
