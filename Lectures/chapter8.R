# bootstrap variance estimation for a median
# Simulate some binary data:
set.seed(1) # repeat yourself for a different random seed
n = 20 
p = .4
x <- rbinom(n=n,size=1,prob=p)
x
T <- median(x)
T
# Using R's default algorithm for the median, the median of 20
# binary observations is 0 if there are 9 or fewer 1's, 0.5 if 
# there are 10 1's and 1 if there are 11 or more 1's. The number
# of 1's in a sample is binomial with success probability p.
# The chance that a sample's median is 0 is
p0 = pbinom(n/2-1,size=n,prob=p)
# the chance that a sample's median is 0.5 is
p0.5 = dbinom(n/2,size=n,prob=p)
# and the chance of 1 is
p1 = 1-pbinom(n/2,size=n,prob=p)
# Variance of T under F:
Tprobs <- c(p0,p0.5,p1)
Tsuppt <-c(0,0.5,1)
EF = sum(Tsuppt*Tprobs)
VF = sum((Tsuppt-EF)^2*Tprobs)
VF

# For the bootstrap, we take the sample as the population. In this
# population, the observed proportion, pobs, is the "true" p.
pobs = mean(x) # proportion in "bootstrap" population
p0boot = pbinom(n/2-1,size=n,prob=pobs)
# the chance that a boostrap sample's median is 0.5 is
p0.5boot = dbinom(n/2,size=n,prob=pobs)
# and the chance of 1 is
p1boot = 1-pbinom(n/2,size=n,prob=pobs)
# Variance of T under F-hat_{T_n}:
Tprobsboot <- c(p0boot,p0.5boot,p1boot)
EFhat = sum(Tsuppt*Tprobsboot)
VFhat = sum((Tsuppt-EFhat)^2*Tprobsboot)
VFhat 
# compare to VF:
VF

# In more complex problems, we would not be able to determine 
# the distribution of T under Fhat, so we might estimate it 
# by sampling.
boot = function(B,x) {
  Tboot <- vector(length=B)
  for(i in 1:B) {
    xstar <- sample(x,size=n,replace=TRUE)
    Tboot[i] <- median(xstar)
  }
  return(Tboot)
}

B = 1000
Tboot = boot(B,x)
table(Tboot)/B
vboot = var(Tboot)
vboot 
# compare to VFhat
VFhat
# already small compared to the difference between VF and VFhat
VF

# bootstrap confidence intervals
# Use the LSAT/GPA data from the text
LSAT<- c(576,635,558,578,666,580,555,661,
         651,605,653,575,545,572,594)
GPA <- c(3.39,3.30,2.81,3.03,3.44,3.07,3.00,3.43,
         3.36,3.13,3.12,2.74,2.76,2.88,2.96) # last entry was typo in text
dat = data.frame(LSAT,GPA)
rm(LSAT,GPA)
T = with(dat,cor(LSAT,GPA))

B = 1000
Tboot = vector(length=B)
for(i in 1:B) {
  inds = sample(1:nrow(dat),replace=TRUE)
  bdat = dat[inds,]
  Tboot[i] = with(bdat,cor(LSAT,GPA))
}
hist(Tboot)
# 1. Normal interval
c(T-sd(Tboot)*qnorm(0.975),T+sd(Tboot)*qnorm(0.975))
# or could truncate upper limit to 1
# 2. Pivot interval
qH = quantile(Tboot-T,probs=c(.975,.025))
T-qH
# 3. Percentile interval
quantile(Tboot,probs=c(.025,.975))

# Bootstrap failure
n=1000 # large n
B = 1000
xnstar = vector(length=B)
x = runif(n)
xn = max(x)
for(j in 1:B) {
  xnstar[j] = max(sample(x,size=n,replace=TRUE))
}
xx = sort(xnstar)-xn
plot(xx,(1:B)/B,type="l") 
lines(xx,(1+xx)^n,col="red")
# repeat for n = 10000