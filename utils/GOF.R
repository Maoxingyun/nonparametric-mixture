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
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result$p.value)
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = result$estimate)
      }
    }
  }
  return(dist_hash)
}


## KL Divergence with degree of membership
GOF_KL_soft = function(X, all_dist, degree_of_membership, threshold){
  
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
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
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "plnorm", mle$estimate[1], mle$estimate[2])
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pexp", mle$estimate[1])
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pgamma", mle$estimate[1], mle$estimate[2])
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result$p.value)
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}



## Kolmogorov-Smirnov Test with degree of membership
GOF_KS_soft = function(X, all_dist, degree_of_membership, threshold){
  
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
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
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "plnorm", mle$estimate[1], mle$estimate[2])
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        result = ks.test(x = data, "pexp", mle$estimate[1])
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
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



## Cram¨¦r¨Cvon Mises criterion
library('CDFt')
GOF_CvM = function(X, all_dist, cluster_lab){
  num_of_samples = 1000
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
        y <- rnorm(num_of_samples, mle$estimate[1], mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        y <- rlnorm(num_of_samples, mle$estimate[1], mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        y <- rexp(num_of_samples, mle$estimate[1])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        y <- rgamma(num_of_samples, shape = mle$estimate[1], scale = 1 / mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result.pvalue)
      
      if (result.pvalue > pvalue) {
        pvalue = result.pvalue
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}


## Cram¨¦r¨Cvon Mises criterion with degree of membership
GOF_CvM_soft = function(X, all_dist, degree_of_membership, threshold){
  num_of_samples = 1000
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        y <- rnorm(num_of_samples, mle$estimate[1], mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        y <- rlnorm(num_of_samples, mle$estimate[1], mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        y <- rexp(num_of_samples, mle$estimate[1])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        y <- rgamma(num_of_samples, shape = mle$estimate[1], scale = 1 / mle$estimate[2])
        result <- CramerVonMisesTwoSamples(data,y)
        result.pvalue = 1/6*exp(-result)
        
      }
      
      if (result.pvalue > pvalue) {
        pvalue = result.pvalue
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}


## Chi square test
library('zoo')
GOF_Chi = function(X, all_dist, cluster_lab){
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
        p1 = hist(x = data, breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pnorm(p1$breaks, mle$estimate[1], mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- plnorm(p1$breaks, mle$estimate[1], mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pexp(p1$breaks, mle$estimate[1])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pgamma(p1$breaks, shape = mle$estimate[1], scale = 1 / mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(a$p.value)
      
      if (a$p.value > pvalue) {
        pvalue = a$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}



## Chi square test with degree of membership
GOF_Chi_soft = function(X, all_dist, degree_of_membership, threshold){
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        p1 = hist(x = data, breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pnorm(p1$breaks, mle$estimate[1], mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- plnorm(p1$breaks, mle$estimate[1], mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pexp(p1$breaks, mle$estimate[1])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        p1 = hist(x = data[which(data>0)], breaks=50, include.lowest=FALSE, right=FALSE, plot=FALSE)
        breaks_cdf <- pgamma(p1$breaks, shape = mle$estimate[1], scale = 1 / mle$estimate[2])
        null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
        a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
        
      }
      
      if (a$p.value > pvalue) {
        pvalue = a$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}



## goftest
library('goftest')
GOF_goftest = function(X, all_dist, cluster_lab, testStr){
  ## testStr = "CvM" / "AD"
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
        if (testStr == "CvM") {
          # result = cvm.test(data, "pnorm", mle$estimate[1], mle$estimate[2], estimated=TRUE)
          # When estimated = TRUE, the test statistics is not stable
          result = cvm.test(data, "pnorm", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "pnorm", mle$estimate[1], mle$estimate[2])
        }

      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "plnorm", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "plnorm", mle$estimate[1], mle$estimate[2])
        }
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "pexp", mle$estimate[1])
        } else {
          result = ad.test(data, "pexp", mle$estimate[1])
        }
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "pgamma", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "pgamma", mle$estimate[1], mle$estimate[2])
        }
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result$p.value)
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}

## goftest with degree of membership
GOF_goftest_soft = function(X, all_dist, degree_of_membership, testStr, threshold){
  ## testStr = "CvM" / "AD"
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          # result = cvm.test(data, "pnorm", mle$estimate[1], mle$estimate[2], estimated=TRUE)
          # When estimated = TRUE, the test statistics is not stable
          result = cvm.test(data, "pnorm", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "pnorm", mle$estimate[1], mle$estimate[2])
        }
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data[which(data>0)], "dlnorm", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "plnorm", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "plnorm", mle$estimate[1], mle$estimate[2])
        }
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data[which(data>0)], "dexp", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "pexp", mle$estimate[1])
        } else {
          result = ad.test(data, "pexp", mle$estimate[1])
        }
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data[which(data>0)], "dgamma", extend = TRUE, relax = TRUE)
        if (testStr == "CvM") {
          result = cvm.test(data, "pgamma", mle$estimate[1], mle$estimate[2])
        } else {
          result = ad.test(data, "pgamma", mle$estimate[1], mle$estimate[2])
        }
        
      }
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}


## goft
library('goft')
GOF_goft = function(X, all_dist, cluster_lab){
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
        data = sample(data) # random shuffle
        if (length(data) > 400) {
          result = normal_test(x = data[1:400]) # norm_test can only process sample size <= 400
        } else {
          result = normal_test(x = data)
        }
        # result = shapiro.test(x = data) # not accurate
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        result = lnorm_test(x = data)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        result = exp_test(x = data, method = "transf") # method can also be "ratio"
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        result = gamma_test(x = data)
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result$p.value)
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}


## gof with degree of membership
GOF_goft_soft = function(X, all_dist, degree_of_membership, threshold){
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        data = sample(data) # random shuffle
        if (length(data) > 400) {
          result = normal_test(x = data[1:400]) # norm_test can only process sample size <= 400
        } else {
          result = normal_test(x = data)
        }
        # result = shapiro.test(x = data) # not accurate
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        result = lnorm_test(x = data)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        result = exp_test(x = data, method = "transf") # method can also be "ratio"
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        result = gamma_test(x = data)
        
      }
      
      if (result$p.value > pvalue) {
        pvalue = result$p.value
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}



## dbEmpLikeGOF
library('dbEmpLikeGOF')
GOF_dbEmpLikeGOF = function(X, all_dist, cluster_lab){
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
        Fx = pnorm(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
      
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        Fx = plnorm(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        Fx = pexp(data, mle$estimate[1])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        Fx = pgamma(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(a$pvalue)
      
      if (a$pvalue > pvalue) {
        pvalue = a$pvalue
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}




## dbEmpLikeGOF with degree of membership
GOF_dbEmpLikeGOF_soft = function(X, all_dist, degree_of_membership, threshold){
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    pvalue = -1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        Fx = pnorm(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        Fx = plnorm(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        Fx = pexp(data, mle$estimate[1])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        Fx = pgamma(data, mle$estimate[1], mle$estimate[2])
        a <- dbEmpLikeGOF(x = Fx, testcall = "uniform", pvl.Table = FALSE, num.mc = 5000)
        
      }
      
      if (a$pvalue > pvalue) {
        pvalue = a$pvalue
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}





## Helliger distance
library('distrEx')
GOF_distrEx = function(X, all_dist, cluster_lab, testStr){
  ## testStr = "Helliger" / "TV"
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
    mindistance = 1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Norm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Norm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        }
        
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Lnorm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Lnorm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        }
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Exp(mle$estimate[1]), data)
        } else {
          result <- TotalVarDist(Exp(mle$estimate[1]), data)
        }
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Gammad(shape=mle$estimate[1], scale=1/mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Gammad(shape=mle$estimate[1], scale=1/mle$estimate[2]), data)
        }
        
      }
      
      # print(i)
      # print(all_dist[[as.character(j)]])
      # print(result)
      
      if (result < mindistance) {
        mindistance = result
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}




## Helliger distance with degree of membership
GOF_distrEx_soft = function(X, all_dist, degree_of_membership, testStr, threshold){
  ## testStr = "Helliger" / "TV"
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    mindistance = 1
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Norm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Norm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        }
        
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Lnorm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Lnorm(mean=mle$estimate[1], sd=mle$estimate[2]), data)
        }
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Exp(mle$estimate[1]), data)
        } else {
          result <- TotalVarDist(Exp(mle$estimate[1]), data)
        }
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        if (testStr == "Helliger") {
          result <- HellingerDist(Gammad(shape=mle$estimate[1], scale=1/mle$estimate[2]), data)
        } else {
          result <- TotalVarDist(Gammad(shape=mle$estimate[1], scale=1/mle$estimate[2]), data)
        }
        
      }
      
      if (result < mindistance) {
        mindistance = result
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}





## Wasserstein distance
library('transport')
GOF_wasserstein = function(X, all_dist, cluster_lab, p){
  # p is the parameter of wasserstein distance
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
    mindistance = 100
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rnorm(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        

      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rlnorm(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rexp(10*length(data), mle$estimate[1]), p = p, wa = NULL, wb = NULL)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rgamma(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        
      }
      
      print(i)
      print(all_dist[[as.character(j)]])
      print(result)
      
      if (result < mindistance) {
        mindistance = result
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}


## Wasserstein distance with degree of membership
GOF_wasserstein_soft = function(X, all_dist, degree_of_membership, p, threshold){
  # p is the parameter of wasserstein distance
  # threshold is about the degree of membership
  
  cluster_hash <- new.env() ## <cluster_lab, data points>
  for (i in 1:ncol(degree_of_membership)) {
    single_X = numeric(0)
    for (j in 1:nrow(degree_of_membership)) {
      if (degree_of_membership[j,i] > threshold) {
        single_X = c(single_X, X[j])
      }
    }
    cluster_hash[[as.character(i)]] <- single_X
  }
  
  
  dist_hash <- new.env() ## <cluster_lab, distribution> 
  ## distribution is a list, name = "distribution name", parameter = c(parametric vector)
  for(i in 1:length(cluster_hash)){
    data = cluster_hash[[as.character(i)]]
    significant_level = 0.05
    mindistance = 100
    
    ## process data, leave only positive number
    processed_data = data[data>0]
    
    ## fit the distribution
    for (j in 1:length(all_dist)) {
      if (all_dist[[as.character(j)]] == "norm") {
        mle = vs.test(x = data, "dnorm", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rnorm(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        
        
      } else if (all_dist[[as.character(j)]] == "lnorm") {
        mle = vs.test(x = data, "dlnorm", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rlnorm(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        
      } else if (all_dist[[as.character(j)]] == "exp") {
        mle = vs.test(x = data, "dexp", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rexp(10*length(data), mle$estimate[1]), p = p, wa = NULL, wb = NULL)
        
      } else if (all_dist[[as.character(j)]] == "gamma") {
        mle = vs.test(x = data, "dgamma", extend = TRUE, relax = TRUE)
        result <- wasserstein1d(data, rgamma(10*length(data), mle$estimate[1], mle$estimate[2]), p = p, wa = NULL, wb = NULL)
        
      }
      
      if (result < mindistance) {
        mindistance = result
        dist_hash[[as.character(i)]] <- list(name = all_dist[[as.character(j)]], parameter = mle$estimate)
      }
    }
  }
  return(dist_hash)
}




