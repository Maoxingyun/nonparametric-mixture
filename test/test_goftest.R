library('goftest')
x <- rnorm(30, mean=2, sd=1)
# default behaviour: parameters fixed: simple null hypothesis
cvm.test(x, "pnorm", mean=2, sd=1)
ad.test(x, "pnorm", mean=2, sd=1)
# parameters estimated: composite null hypothesis
mu <- mean(x)
sigma <- sd(x)
result1 = cvm.test(x, "pnorm", mean=mu, sd=sigma, estimated=TRUE)
result2 = ad.test(x, "pnorm", mean=mu, sd=sigma, estimated=TRUE)



pAD(0.77, n=5)
pAD(0.77)
pAD(0.77, fast=FALSE)
qAD(0.5, n=5)
qAD(0.5)

pCvM(0.12, n=5)
pCvM(0.12)
qCvM(0.5, n=5)
qCvM(0.5)


data <- rnorm(300, 10,1)
testStr = "CvM"