# bootstrap variance estimation for a median
# Simulate some binary data:
set.seed(1)
n = 20
p = .3
x <- rbinom(n=n,size=1,prob=p)
x
pobs = mean(x) # proportion in "bootstrap" population
T <- median(x)
T
# Using R's default algorithm for the median, the median of 20
# binary observations is 0 if there are 9 or fewer 1's, 0.5 if 
# there are 10 1's and 1 if there are 11 or more 1's. The number
# of 1's in a sample is binomial with success probability p.
# The chance that a sample's median is 0 is
p0 = pbinom(9,size=20,prob=p)
# the chance that a sample's median is 0.5 is
p0.5 = dbinom(10,size=20,prob=p)
# and the chance of 1 is
p1 = 1-pbinom(10,size=20,prob=p)
# Variance of T under F:
Tprobs <- c(p0,p0.5,p1)
Tsuppt <-c(0,0.5,1)
EF = sum(Tsuppt*Tprobs)
VF = sum((Tsuppt-EF)^2*Tprobs)
VF

# For the bootstrap, we take the sample as the population. In this
# population, the observed proportion, pboot, is the "true" p.
p0boot = pbinom(9,size=20,prob=pobs)
# the chance that a boostrap sample's median is 0.5 is
p0.5boot = dbinom(10,size=20,prob=pobs)
# and the chance of 1 is
p1boot = 1-pbinom(10,size=20,prob=pobs)
# Variance of T under F-hat_{T_n}:
Tprobsboot <- c(p0boot,p0.5boot,p1boot)
EFhat = sum(Tsuppt*Tprobsboot)
VFhat = sum((Tsuppt-EFhat)^2*Tprobsboot)
VFhat
VFhat - VF

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
vboot - VFhat
# already small compared to the difference between VF and VFhat
VFhat - VF
