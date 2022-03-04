# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Enhanced mode clustering." arXiv preprint arXiv:1406.1780 (2014).
# Date: 04/07/2015
source("EMC.R")


##### Generating Data
library(freqparcoord)
library('TDA')

# 四个子类
# 每个类100个点
mean1 <- matrix(c(0,5),nrow = 2)
Sigma <- matrix(c(1,0,0,1),2,2)
D = mvrnorm(n = 100, mean1, Sigma)
D = cbind(matrix(1, nr = 100, nc = 1), D)

mean2 <- matrix(c(5,0),nrow = 2)
D = rbind(D, cbind(matrix(2, nr = 100, nc = 1), mvrnorm(n = 100, mean2, Sigma)))

mean3 <- matrix(c(0,-5),nrow = 2)
D = rbind(D, cbind(matrix(3, nr = 100, nc = 1), mvrnorm(n = 100, mean3, Sigma)))

mean4 <- matrix(c(-5,0),nrow = 2)
D = rbind(D, cbind(matrix(4, nr = 100, nc = 1), mvrnorm(n = 100, mean4, Sigma)))

# 可加outliers
head(D)

X0= D[,2:3]
X = scale(X0)
plot(X, xlim = c(-4,4),ylim = c(-4,4), main="Scaled original data")
Y = as.numeric(D[,1])


## (1.1) Nonarametric clustering: Enhanced mode clustering
X_emc = EMC(X)
nonpara_result = X_emc$labels

## (1.2) Nonarametric clustering: "Density Clustering" in library('TDA')
## 可以聚类一维数据
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
plot(X, xlim = c(-4,4),ylim = c(-4,4), main="Results of Nonparametric Cluster", col=rainbow(max(nonpara_result))[nonpara_result])



## (2.1) Goodness of fit (Google: Goodness of fit in R)
## 'gofTest' in R
## library('goft' / 'vsgoftest')
## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/





