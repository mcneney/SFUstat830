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


# Test the null hypothesis that the variance is 1 with
# a Wald test with reference distribution from 
# (i) the asymptotic normality, or 
# (ii) parametric bootstrap.
#-----------------------------------------------------------------#
#-----------------------------------------------------------------#
# 1. Asymptotic normality -- three possible SEs
#-----------------------------------------------------------------#
# 1.1 SE from the observed Fisher information, which is diagonal,
#-----------------------------------------------------------------#
# so the variance is estimated by the lower-right element of
# [-\ddot\l(\hat\theta)]^{-1} = \hat\sigma^4/n.
# The test statistic is (\hat\sigma^2-1)/\sqrt{\hat\sigma^4/n}
teststat = function(x,null=1) {
  n = length(x); ss = var(x)*(n-1)/n
  return((ss-null)/(ss/sqrt(n))) 
}
out = runsim(Nsim,mean=0) # doesn't matter what mean is
hist(out$teststats,nclass=30) # distribution is skewed
hist(out$pvals,nclass=30) 
# p-value distribution doesn't look uniform, test is anti-conservative
mean(out$pvals<=0.05) # size of about 0.19 > 0.05
#-----------------------------------------------------------------#
# 1.2 SE from the Fisher information evaluated at the MLE of mu and
#-----------------------------------------------------------------#
# the hypothesized value of sigma^2.
# Can show that under (\bar{X}_n,\sigma^2) the variance of \hat\sigma^2 is 
# (2 \sigma^6)/n * {2\hat\sigma^2 - \sigma^2}^{-1}
teststat = function(x,null=1) {
  n = length(x); ss = var(x)*(n-1)/n
  se = sqrt(2*null^3/(n*(2*ss-null))) # var can be <0 so sqrt can be NaN
  return((ss-null)/se) 
}
out = runsim(Nsim,mean=0)
hist(out$teststats,nclass=30)
hist(out$pvals,nclass=30) # distribution is non-uniform
mean(out$pvals<=0.05,na.rm=TRUE) # size of about 0.19 > 0.05
# Less bias in the test, but a very strange p-value distribution.
#-----------------------------------------------------------------#
# 1.3 SE by bootstrap. I used parametric bootstrap.
#-----------------------------------------------------------------#
parboot.pval = function(x,Nboot=1000) {
  xbar = mean(x); n = length(x); ss = var(x)*(n-1)/n
  tt = ss - 1
  tstar = vector(length=Nboot)
  for(j in 1:Nboot) {
    xstar = simdat(n,xbar,sqrt(ss))
    tstar[j] = var(xstar)*(n-1)/n - ss # in bootstrap world, \hat\sigma^2 is truth
  }
  se = sd(tstar)
  Z = tt/se
  pp = 2*pnorm(-abs(Z))
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
Nsim = 1000 # Can't afford 10000 simulations
system.time({out = runsim.parboot(Nsim,mean=0)})# about 40 sec
hist(out$pvals,nclass=30) #  more uniform than Fisher-info-based approaches
mean(out$pvals <= 0.05) # but still anti-conservative
#-----------------------------------------------------------------#
#-----------------------------------------------------------------#
# (2) parametric bootstrap
#-----------------------------------------------------------------#
# Generate a bootstrap distribution for the Wald test statistic by sampling from 
# N(\hat\mu,\hat\sigma^2). Use the Wald test statistic based on the
# observed Fisher information, since this one didn't have problems with
# negative variances.
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
# Revert to the observed Fisher info statistic, which didn't have problems with 
# taking square-root of negative variances.
teststat = function(x,null=1) {
  n = length(x); ss = var(x)*(n-1)/n
  return((ss-null)/(ss/sqrt(n))) 
}
Nsim = 1000
system.time({out = runsim.parboot(Nsim,mean=0)})# about 45 sec
hist(out$pvals,nclass=30) # much more uniform than the ML approaches
mean(out$pvals <= 0.05) # about right
#-----------------------------------------------------------------#
#-----------------------------------------------------------------#
# Other things to try: 
# - What happens to ML as the sample size increases?
# - What if the data-generating distribution was not normal?
#   Then how would the parametric bootstrap work? 
#   Would nonparametric bootstrap work?