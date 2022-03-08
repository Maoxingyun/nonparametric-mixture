library("TDA")# 四个子类
library('vsgoftest')
library(MASS)
source("EMC_mod.R")
# 每个类100个点

D = rnorm(n = 500, -5,1)
D = cbind(matrix(1, nr = 500, nc = 1), D)
D = rbind(D, cbind(matrix(2, nr = 300, nc = 1), rgamma(n = 300, 24, 4)))
D = rbind(D, cbind(matrix(3, nr = 400, nc = 1), rexp(n = 400, 2)))

X0 = D[,2:2]
X = scale(X0) ## scale only for EMC nonparametric clustering

yy = density(X0)
plot(yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
points(X0, y=rep(0,nrow(X)), pch=1)


## (1) Nonparametric clustering
## (1.1) Enhanced mode clustering
X_emc = EMC(X)
cluster_lab = X_emc$labels


# ## (1.2) "Density Clustering" in library('TDA')
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




cluster_hash <- new.env() ## <cluster_lab, data points>
for(i in 1:max(cluster_lab)){
  single_X = X0[which(cluster_lab==i)]
  cluster_hash[[as.character(i)]] <- single_X
}


dist_hash <- new.env() ## <cluster_lab, distribution> 
## distribution is a list, name = "distribution name", parameter = c(parametric vector)
for(i in 1:max(cluster_lab)){
  data = cluster_hash[[as.character(i)]]
  significant_level = 0.05
  pvalue = 0
  
  ## process data, leave only positive number
  processed_data = data[data>0]
  
  ## fit the distribution
  for (j in 1:total_dist) {
    if (all_dist[[as.character(j)]] != "norm") {
      if (length(processed_data) < 0.5*length(data)) {
        next
      }
      result = vs.test(x = processed_data, densfun = paste("d",all_dist[[as.character(j)]],sep=""))
    } else {
      result = vs.test(x = data, densfun = paste("d",all_dist[[as.character(j)]],sep=""))
    }
    
    if (result$p.value > pvalue) {
      pvalue = result$p.value
      dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = result$estimate)
    }
  }
}



## (2.1) ## xxx.test see https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/



## (2.2) library('goft' / 'vsgoftest')
## (2.3) 'gofTest' in R

## (2.4) Methods in "All of statistics"




