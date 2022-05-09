num_of_samples = 1000
x <- rgamma(num_of_samples, shape = 10, scale = 3)
x <- x + rnorm(length(x), mean=0, sd = .1)

p1 <- hist(x,breaks=50, include.lowest=FALSE, right=FALSE)

# Chi Square test
library('zoo')
breaks_cdf <- pgamma(p1$breaks, shape=10, scale=3)
null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)

# Cram¨¦r¨Cvon Mises criterion
library('CDFt')
num_of_samples = 100000
y <- rgamma(num_of_samples, shape = 10, scale = 3)
res <- CramerVonMisesTwoSamples(x,y)
p.value = 1/6*exp(-res)

# Kolmogorov¨CSmirnov test
result = ks.test(x, "pgamma", 10, 1/3)
