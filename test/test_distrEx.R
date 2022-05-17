library('distrEx')

test = Norm(mean = 50, sd = sqrt(25))

x <- rnorm(400, 10, 1)
res_x <- HellingerDist(Norm(mean = 10, sd = 1), x)

y <- rlnorm(400, 0, 1)
res_y <- HellingerDist(Lnorm(0, 1), y)

p <- rexp(400, 3)
res_p <- HellingerDist(Exp(3), p)


q <- rgamma(400, 24, 4)
res_q <- HellingerDist(Gammad(scale=0.25, shape=24), q)


