library("TDA")# 四个子类
library(MASS)
source("EMC_mod.R")
# 每个类100个点
D = rnorm(n = 100, 6, 0.25)
D = cbind(matrix(1, nr = 100, nc = 1), D)
D = rbind(D, cbind(matrix(2, nr = 100, nc = 1), rnorm(n = 100, 2, 0.25)))
D = rbind(D, cbind(matrix(3, nr = 100, nc = 1), mvrnorm(n = 100, -2, 0.25)))
D = rbind(D, cbind(matrix(4, nr = 100, nc = 1), mvrnorm(n = 100, -6, 0.25)))

X0 = D[,2:2]
X = scale(X0)


yy = density(X)
plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
points(X, y=rep(0,nrow(X)), pch=1)



## (1.1) Nonarametric clustering: Enhanced mode clustering
X_emc = EMC(X)
nonpara_result = X_emc$labels


## (1.2) Nonarametric clustering: "Density Clustering" in library('TDA')
Tree <- clusterTree(X, k = 100, density = "knn",printProgress = FALSE)
TreeKDE <- clusterTree(X, k = 100, h = 0.3, density = "kde",printProgress = FALSE)
plot(Tree, type = "lambda", main = "lambda Tree (knn)")
plot(Tree, type = "kappa", main = "kappa Tree (knn)")
plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")
plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")

nonpara_result <- vector(length = 400)
for (i in 1:length(TreeKDE$DataPoints)) {
  for (j in TreeKDE$DataPoints[i][[1]]) {
    nonpara_result[j] = i
  }
}

## Visualization
plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
points(X, y=rep(0,nrow(X)), pch=1, col=rainbow(max(nonpara_result))[nonpara_result])


## (2.1) Goodness of fit (Google: Goodness of fit in R)
## 'gofTest' in R
## library('goft' / 'vsgoftest')
## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/

