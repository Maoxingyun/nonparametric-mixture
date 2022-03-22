# library("TDA")
# library('vsgoftest')
library(MASS)
source("Clustering.R")
source("GOF.R")
source("Refinement.R")

D = rnorm(n = 500, -5,1)
D = cbind(matrix(1, nr = 500, nc = 1), D)
D = rbind(D, cbind(matrix(2, nr = 300, nc = 1), rgamma(n = 300, 24, 4)))
D = rbind(D, cbind(matrix(3, nr = 400, nc = 1), rexp(n = 400, 2)))

X0 = D[,2:2]

yy = density(X0)
colscale = c("black", "red")
plot(yy, type="l", ylab="Density", xlab="x", ylim=c(0,1.5*max(yy$y)), las=1, lwd=2, main="Nonparametric density estimation")

xx = seq(-10, 15, length=200)
real_density = (5/12)*dnorm(xx, -5, 1) + (3/12)*dgamma(xx, 24, 4) + (4/12)*dexp(xx, 2)
lines(xx, real_density, col=colscale[2], lty=2, lwd=2)

points(X0, y=rep(0,length(X0)), pch=1)
legend(7, 0.10, c("KDE","True distribution"), col=colscale[c(1,2)], lty=c(1,2), lwd=2, bty="n")

max_iteration = 1
for (i in 1:max_iteration) {
  X = scale(X0) ## scale only for EMC nonparametric clustering
  ## (1) Nonparametric clustering
  cluster_lab = Clustering(X) ## EMC method
  
  # TODO
  ## (1.1) "Density Clustering" in library('TDA')
  # Tree <- clusterTree(X, k = 100, density = "knn",printProgress = FALSE)
  # TreeKDE <- clusterTree(X, k = 100, h = 0.3, density = "kde",printProgress = FALSE)
  # plot(Tree, type = "lambda", main = "lambda Tree (knn)")
  # plot(Tree, type = "kappa", main = "kappa Tree (knn)")
  # plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")
  # plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")
  # 
  # cluster_lab <- vector(length = nrow(X))
  # for (i in 1:length(TreeKDE$DataPoints)) {
  #   for (j in TreeKDE$DataPoints[i][[1]]) {
  #     cluster_lab[j] = i
  #   }
  # }
  ## (1.2) log-concave
  
  ## Visualization
  plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
  points(X0, y=rep(0,nrow(X)), pch=1, col=rainbow(max(cluster_lab))[cluster_lab])
  
  
  
  ## (2) Goodness of fit (Google: Goodness of fit in R)
  all_dist <- new.env()
  all_dist[["1"]] <- "norm"
  all_dist[["2"]] <- "lnorm"
  all_dist[["3"]] <- "exp"
  all_dist[["4"]] <- "gamma"
  all_dist[["5"]] <- "f"
  total_dist = 5
  
  
  dist_hash = GOF(X0, all_dist, cluster_lab) ## vsgoftest method
  
  
  ## TODO
  ## (2.1) ## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
  ## (2.2) library('goft')
  ## (2.3) 'gofTest' in R
  ## (2.4) Methods in "All of statistics"
  
  ## (3) Data refinement
  X0 = Refinement(X0, dist_hash, cluster_lab)
  
  ## Visualization
  yy = density(X0)
  plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
}



