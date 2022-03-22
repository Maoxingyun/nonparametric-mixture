#' Data Refinement
#' 
#' @param X Input data matrix.
#' @param dist_hash result of goodness of fit
#' @param cluster_lab result of nonparametric clustering
#' @return Y: refined data based on dist_hash
#'  
#' @export
#' 
#' 
#' 
ComputeDensity = function(x, dist, parameter){
  if(dist == "norm"){
    return (dnorm(x, parameter[1], parameter[2]))
  } else if(dist=="lnorm"){
    return (dlnorm(x, parameter[1], parameter[2]))
  } else if(dist=="exp"){
    return (dexp(x, parameter[1]))
  } else if(dist=="gamma"){
    return (dgamma(x, parameter[1], parameter[2]))
  } else if(dist=="f"){
    return (df(x, parameter[1], parameter[2]))
  }
  return (0)
}

Refinement = function(X, dist_hash, cluster_lab){
  Y = numeric(0)
  threshold = 0.01
  
  for(i in 1:max(cluster_lab)){
    single_X = X[which(cluster_lab==i)]
    ## process single_X with dist_hash[[as.character(i)]] here
    single_X <- single_X[ ComputeDensity(single_X, dist_hash[[as.character(i)]]$name, dist_hash[[as.character(i)]]$parameter) >= threshold ]
    Y = c(Y, single_X)
  }
  return(Y)
}