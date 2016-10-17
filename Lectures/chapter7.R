dat = scan("http://www.stat.cmu.edu/~larry/all-of-statistics/=data/nerve.dat")

ee = ecdf(dat)
plot(ee)
