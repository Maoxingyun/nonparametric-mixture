## (2.1) ## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
## p-value > significant level, 0.01 <= significant level <= 0.1
## Chi Square test
library('zoo')

num_of_samples = 1000
x <- rgamma(num_of_samples, shape = 10, scale = 3)
x <- x + rnorm(length(x), mean=0, sd = .1)

p1 <- hist(x,breaks=50, include.lowest=FALSE, right=FALSE) ## Primary distribution


breaks_cdf <- pgamma(p1$breaks, shape=10, scale=3) ## Reference distribution
null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)

## Cram¨¦r¨Cvon Mises criterion
library('CDFt')
y <- rgamma(num_of_samples, shape = 10, scale = 3)
res <- CramerVonMisesTwoSamples(x,y)
p_value = 1/6*exp(-res)


## Kolmogorov¨CSmirnov test
y <- rgamma(num_of_samples, shape = 10, scale = 3)
result = ks.test(x, y)

## (2.2) 'vsgoftest' in R
library('vsgoftest')

## entropy.estimate
set.seed(2) #set seed of PRNG
samp <- rnorm(n = 100, mean = 0, sd = 1) #sampling from normal distribution
entropy.estimate(x = samp, window = 8) #estimating entropy with window = 8
log(2*pi*exp(1))/2 #the exact value of entropy

## vs.test
set.seed(5)
samp <- rnorm(50,2,3)
vs.test(x = samp, densfun = 'dlaplace')

set.seed(4)
vs.test(x = samp, densfun = 'dnorm')

set.seed(26)
vs.test(x = samp, densfun = 'dnorm', param = c(2,3))

set.seed(1)
samp <- rweibull(200, shape = 1.05, scale = 1)
set.seed(2)
vs.test(samp, densfun = 'dexp', simulate.p.value = TRUE, B = 10000)

set.seed(63)
vs.test(samp, densfun = 'dexp', delta = 5/30)

set.seed(8)
samp <- rexp(30, rate = 3)
vs.test(x = samp, densfun = "dlnorm")
vs.test(x = samp, densfun = "dlnorm", extend = TRUE)

samp <- c(samp, rep(4,3)) #add ties in the previous sample
vs.test(x = samp, densfun = "dexp")
vs.test(x = samp, densfun = "dexp", extend = TRUE)





