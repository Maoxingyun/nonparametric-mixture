#' Nonparametric clustering
#' 
#' @param X Input data matrix.
#' @return cluster_lab: result of nonparametric clustering
#'  
#' @export
source('EMC_mod.R')
Clustering = function(X){
  X_emc = EMC(X)
  cluster_lab = X_emc$labels
  return(cluster_lab)
}