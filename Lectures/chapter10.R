rm(list=ls())

# Chapter 10, hypothesis testing and p-values
# Under H_0, the distribution of the p-value is uniform on [0,1]

n = 50
simdat = function(n,mean,sd=1) {
  return(rnorm(n,mean=mean,sd=sd))
}
teststat = function(x,null=0,sd=1) {
  xbar = mean(x); n = length(x)
  return((xbar-null)/(sd/sqrt(n))) 
}
pval = function(t) {
  return(2*pnorm(-abs(t)))
}

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
# The test statistic is now (\hat\sigma^2-1)/\sqrt{\hat\sigma^4/n}
teststat = function(x,null=1) {
  n = length(x); ss = var(x)*(n-1)/n
  return((ss-null)/(ss/sqrt(n))) # testing H_0: sigma^2=1, vs H_1: sigma^2 not= 1
}
out = runsim(Nsim,mean=0) # doesn't matter what mean is
hist(out$teststats,nclass=30) # distribution is skewed
hist(out$pvals,nclass=30) 
# p-value distribution doesn't look uniform, test is anti-conservative
mean(out$pvals<=0.05) # size of about 0.19 > 0.05
#-----------------------------------------------------------------#
# (2) parametric bootstrap
# We generate a bootstrap distribution for the test statistic by sampling from 
# N(\hat\mu,\hat\sigma^2).
parboot.pval = function(x,Nboot=1000) {
  xbar = mean(x); n = length(x); ss = var(x)*(n-1)/n
  tt = teststat(x,null=1)
  tstar = vector(length=Nboot)
  for(j in 1:Nboot) {
    xstar = simdat(n,xbar,sqrt(ss))
    tstar[j] = teststat(xstar,null=ss) # in bootstrap world, \hat\sigma^2 is truth
  }
  pp = mean((tstar>=abs(tt)) | (tstar<= (-abs(tt))))
  return(pp)
}
runsim.parboot = function(Nsim,mean) {
  pp = vector(length=Nsim)
  for(i in 1:Nsim) {
    dat = simdat(n,mean=mean)
    pp[i] = parboot.pval(dat)
  }
  return(list(pvals = pp))
}
Nsim = 1000
system.time({out = runsim.parboot(Nsim,mean=0)})# about 45 sec
hist(out$pvals,nclass=30)
mean(out$pvals <= 0.05) # about right
