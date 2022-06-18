#' Determine degree of membership from result from nonparametric clustering
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




Membership_OKM = function(X, result){
  Y = matrix(0, nrow = length(X), ncol = length(result$modes))
  
  for (i in 1:length(X)) {
    phi = 0
    cnt = 0 # the number of elements in A_i  
    while (1) {
      if (cnt == length(result$modes)) {
        break
      }
      
      # m_star = argmin_m |xi-m|^2, where m is not in A_i
      m_star = -1
      dist = Inf
      for(k in 1:length(result$modes)) {
        if (Y[i,k] == 0 && (m_star == -1 || abs(X[i]-result$modes[k]) < dist)) {
          dist = abs(X[i]-result$modes[k])
          m_star = k
        }
      }
      Y[i, m_star] = 1
      
      # compute newPhi with A_i+{m_star}
      newPhi = (cnt*phi+result$modes[m_star]) / (cnt+1)
  
      if (cnt == 0 || abs(X[i]-phi) > abs(X[i]-newPhi)) {
        phi = newPhi
        cnt = cnt+1
      } else {
        Y[i, m_star] = 0
        break
      }
    }
  }
  
  return(Y)
}


