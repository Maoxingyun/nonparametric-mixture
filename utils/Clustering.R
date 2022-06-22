#' Nonparametric clustering
#' 
#' @param X Input data matrix.
#' @return cluster_lab: result of nonparametric clustering
#'  
#' @export

source('EMC_mod.R')
Clustering_EMC = function(X0){
  X = scale(X0) ## scale only for EMC nonparametric clustering
  X_emc = EMC(X)
  result <- list(name = "Nonparametric_result", cluster_lab = X_emc$labels, modes = X_emc$modes)
  result$modes = result$modes * sd(X0)+mean(X0)
  return(result)
}

library('reticulate')
use_python("D:\\Anaconda\\python.exe")
source_python("awc.py")
Clustering_AWC = function(X){
  AWC_object <- AWC(speed=1., n_neigh=200, n_outliers = 20)
  np <- import("numpy", convert = FALSE)
  l = 25
  AWC_object$awc(l, np$array(X))
  clusters <- AWC_object$get_clusters()
  cluster_lab <- AWC_object$get_labels()
  cluster_lab = cluster_lab+1
  return(cluster_lab)
}