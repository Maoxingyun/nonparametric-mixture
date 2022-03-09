#' Data Refinement
#' 
#' @param X Input data matrix.
#' @param dist_hash result of goodness of fit
#' @param cluster_lab result of nonparametric clustering
#' @return Y: refined data based on dist_hash
#'  
#' @export
Refinement = function(X, dist_hash, cluster_lab){
  Y = numeric(0)
  for(i in 1:max(cluster_lab)){
    single_X = X[which(cluster_lab==i)]
    ## TODO: process single_X with dist_hash[[as.character(i)]] here

    Y = c(Y, single_X)
  }
  return(Y)
}