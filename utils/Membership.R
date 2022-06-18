#' Data Refinement
#' 
#' @param X Input data matrix.
#' @param result result of nonparametric clustering (S3 object)
#'        @include cluster_lab: partition of nonparametric clustering
#'        @include modes: modes of clusters
#' @return Y: degree of membership matrix (dimension: [length(X) * modes(X)])
#'  
#' @export
#' 
#' 
#' 
#' 
Membership_fuzzy = function(X, result, weight){
  Y = matrix(0, nrow = length(X), ncol = length(result$modes))
  
  for (i in 1:length(X)) {
    # sum_k 1/|xi-mk|^(1/(m-1))
    sum = 0
    for(k in 1:length(result$modes)) {
      sum = sum+1/abs(X[i]-result$modes[k])^(1/(weight-1))
    }
    
    # |xi-mj|^(1/(m-1))
    for(j in 1:length(result$modes)) {
      Y[i,j] = (abs(X[i]-result$modes[j])^(1/(weight-1)) * sum)^(-1)
    }
    
  }

  return(Y)
}