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

## KL Divergence
GOF_KL = function(X, all_dist, cluster_lab){
  
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
    for (j in 1:length(all_dist)) {
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




## Kolmogorov-Smirnov Test
GOF_KS = function(X, all_dist, cluster_lab){
  
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
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pnorm", mle$estimate[1], mle$estimate[2])
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "plnorm", mle$estimate[1], mle$estimate[2])
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pexp", mle$estimate[1])
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pgamma", mle$estimate[1], mle$estimate[2])
      }
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}