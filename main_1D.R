library("TDA")
library(fossil)
library(MASS)
source("Clustering.R")
source("GOF.R")
source("Refinement.R")

# ========== Generate data ========== #
# Sub-distributions
D = rnorm(n = 500, -5,1)
D = cbind(matrix(1, nr = 500, nc = 1), D)
D = rbind(D, cbind(matrix(2, nr = 300, nc = 1), rgamma(n = 300, 24, 4)))
D = rbind(D, cbind(matrix(3, nr = 400, nc = 1), rexp(n = 400, 2)))
# Add Noise
D = rbind(D, cbind(matrix(-1, nr = 300, nc = 1), runif(n = 300, min = -10, max = 10))) # uniform noise
## TODO: Gaussian noise

X0 = D[,2:2]
labels = D[,1:1]
print(paste("Raw data size: ", length(X0)))

# ========== Visualize raw data ========== # 

yy = density(X0)
colscale = c("black", "red")
plot(yy, type="l", ylab="Density", xlab="x", ylim=c(0,1.5*max(yy$y)), las=1, lwd=2, main="Raw data")

xx = seq(-10, 15, length=200)
real_density = (5/12)*dnorm(xx, -5, 1) + (3/12)*dgamma(xx, 24, 4) + (4/12)*dexp(xx, 2)
lines(xx, real_density, col=colscale[2], lty=2, lwd=2)

points(X0, y=rep(0,length(X0)), pch=1)
legend(7, 0.10, c("KDE","True distribution"), col=colscale[c(1,2)], lty=c(1,2), lwd=2, bty="n")


# ========== Clustering algorithm ========== # 

max_iteration = 3
for (i in 1:max_iteration) {
  print(paste("Iteration", i))
  # ========== step 1: Nonparametric clustering ========== #
  print("Step 1...")
  X = scale(X0) ## scale only for EMC nonparametric clustering
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
  plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2, main=paste("Iteration", i, ",Step 1"))
  points(X0, y=rep(0,nrow(X)), pch=1, col=rainbow(max(cluster_lab))[cluster_lab])
  legend(7, 0.10, c("KDE"), col=colscale[c(1)], lty=c(1), lwd=2, bty="n")
  
  # Evaluate the result
  print(paste("rand index:", rand.index(labels[which(labels!=-1)], cluster_lab[which(labels!=-1)])))
  print(paste("adjusted rand index:", adj.rand.index(labels[which(labels!=-1)], cluster_lab[which(labels!=-1)])))
  
  
  # ========== step 2: Goodness of fit ========== # 
  print("Step 2...")
  all_dist <- new.env()
  all_dist[["1"]] <- "norm"
  all_dist[["2"]] <- "lnorm"
  all_dist[["3"]] <- "exp"
  all_dist[["4"]] <- "gamma"
  all_dist[["5"]] <- "f"
  total_dist = 5
  
  
  dist_hash = GOF(X0, all_dist, cluster_lab) ## vsgoftest method
  
  ## Print the result
  for (j in 1:length(dist_hash)) {
    print(paste("Distribution", j,":",dist_hash[[as.character(j)]]$name))
    print(dist_hash[[as.character(j)]]$parameter)
  }
  
  
  ## TODO
  ## (2.1) ## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
  ## (2.2) library('goft')
  ## (2.3) 'gofTest' in R
  ## (2.4) Methods in "All of statistics"
  
  # ========== step 3: Data refinement ========== #
  print("Step 3...")
  refinement_result = Refinement(X0, labels, dist_hash, cluster_lab)
  X0 = refinement_result$data
  labels = refinement_result$labels
  print(paste("Data size after Refinement: ", length(X0)))

  ## Visualization
  yy = density(X0)
  plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2, main=paste("Iteration", i, ",Step 3"))
  points(X0, y=rep(0,length(X0)), pch=1)
  legend(7, 0.10, c("KDE"), col=colscale[c(1)], lty=c(1), lwd=2, bty="n")
  print("")
}



