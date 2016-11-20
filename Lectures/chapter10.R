rm(list=ls())

# Chapter 10, hypothesis testing and p-values
# Under H_0, the distribution of the p-value is uniform on [0,1]

n = 50
simdat = function(n,mean,sd=1) {
  return(rnorm(n,mean=mean,sd=sd))
}

teststat = function(x,sd=1) {
  xbar = mean(x); n = length(x)
  return(xbar/(sd/sqrt(n))) # testing H_0: mean=0, vs H_1: mean not= 0
}
pval = function(t) {
  return(2*pnorm(-abs(t)))
}

Nsim = 10000
runsim = function(Nsim,mean) {
  tt = pp = vector(length=Nsim)
  for(i in 1:Nsim) {
    dat = simdat(n,mean=mean)
    tt[i] = teststat(dat)
    pp[i] = pval(tt[i])
  }
  return(list(teststats = tt,pvals = pp))
}

Nsim = 10000
out = runsim(Nsim,mean=0)
hist(out$teststats,nclass=30)
hist(out$pvals,nclass=30)

# Change the mean in the simulations
out = runsim(Nsim,mean=.1); hist(out$pvals,nclass=30)
out = runsim(Nsim,mean=.2); hist(out$pvals,nclass=30)
out = runsim(Nsim,mean=.3); hist(out$pvals,nclass=30)
out = runsim(Nsim,mean=.4); hist(out$pvals,nclass=30)
out = runsim(Nsim,mean=.5); hist(out$pvals,nclass=30)


# Test the null hypothesis that the variance is 1 by (i) ML,
# (ii) parametric bootstrap, and (iii) nonparametric bootstrap.
# 1. ML
# We saw in class that the MLE \hat\sigma^2 of the variance 
# \sigma^2 is approximately normal with mean \sigma^2 and
# variance that is estimated by the lower-right element of
# [-\ddot\l(\hat\theta)]^{-1} = \hat\sigma^4/n (WHY?).


