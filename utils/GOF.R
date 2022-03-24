#' Goodness of fit
#' 
#' @param X Input data matrix.
#' @param all_dist The available distributions.
#' @param cluster_lab Result of nonparametric clustering.
#' @return dist_hash a hashmap<cluster_lab, distribution>; where distribution a list consisting of
#'  \item{name}{distribution name}
#'  \item{parameter}{a parametric vector}
#'  
#' @export
library('vsgoftest')
GOF = function(X, all_dist, cluster_lab){
  
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for(i in 1:max(cluster_lab)){
    single_X = X[which(cluster_lab==i)]
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:max(cluster_lab)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
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
  return(dist_hash)
}