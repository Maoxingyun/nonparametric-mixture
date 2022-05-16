library('goft')
data("strength")
comp_strength <- strength$cstrength # compressive strength
res = lnorm_test(comp_strength)


library(fitdistrplus)
data("danishuni")
loss <- danishuni$Loss # losses
logloss <- log(loss) # logarithm of losses
logloss <- sort(logloss[logloss > 0]) # only positive observations are kept
gamma_test(logloss)
gam.fit <- gamma_fit(logloss); gam.fit # fitting the gamma distribution
