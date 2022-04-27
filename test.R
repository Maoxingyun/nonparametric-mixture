setwd("D:/HW/¿ÆÑÐ/Code/nonparametric-mixture/utils")
library('reticulate')
use_python("D:\\Anaconda\\python.exe")
source_python("awc.py")
AWC_object <- AWC(speed=1., n_neigh=200)
np <- import("numpy", convert = FALSE)


# ========== Generate data ========== #
# Sub-distributions
D = rnorm(n = 400, 10,1)
D = cbind(matrix(1, nr = 400, nc = 1), D)
D = rbind(D, cbind(matrix(2, nr = 400, nc = 1), rgamma(n = 400, 24, 4)))
D = rbind(D, cbind(matrix(3, nr = 400, nc = 1), rexp(n = 400, 2)))
X0 = D[,2:2]

l = 25
AWC_object$awc(l, np$array(X0))
clusters <- AWC_object$get_clusters()
labels <- AWC_object$get_labels()

