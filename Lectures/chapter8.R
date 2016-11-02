# bootstrap 
set.seed(1)
n = 20
p = .3


x <- rbinom(n=n,size=1,prob=p)
pobs = mean(x)
x
T <- median(x)
T
# Using R's default algorithm for the median, the median of 20
# binary observations is 0 if there are 9 or fewer 1's, 0.5 if 
# there are 10 1's and 1 if there are 11 or more 1's. The number
# of 1's in a boostrap sample from x is binomial with success 
# probability pobs.
# The chance that a boostrap sample's median is 0 is
p0 = pbinom(9,size=20,prob=pobs)
# the chance that a boostrap sample's median is 0.5 is
p0.5 = dbinom(10,size=20,prob=pobs)
# and the chance of 1 is
p1 = 1-pbinom(10,size=20,prob=pobs)

# Variance of T under F-hat_{T_n}:
Tprobs <- c(p0,p0.5,p1)
Tsuppt <-c(0,0.5,1)
EFhat = sum(Tsuppt*Tprobs)
VFhat = sum((Tsuppt-EFhat)^2*Tprobs)
VFhat

# Compare to true variance estimated by simulation from the true model.
Nsim <- 1e5
Tsim <- vector(length=Nsim)
for(i in 1:Nsim) {
  xsim <- rbinom(n=n,size=1,prob=p)
  Tsim[i] <- median(xsim)
}
VFsim = var(Tsim)
VFsim


# Use random sampling to appoximate this distribution.
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
table(Tboot)/B - Fhat.Tn.p
var(Tboot) - VFhat

B = 10000
table(boot(B,x))/B - Fhat.Tn.p

B = 100000
table(boot(B,x))/B - Fhat.Tn.p
