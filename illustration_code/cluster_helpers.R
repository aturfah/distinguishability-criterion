## Functions to implement the different clustering procedures
## Date: 4/20/2024
library(mclust)

source("code/Pmc_implementation.R")

#' Get the hierarchical clustering tree for given data
#' 
#' @param X Vector or matrix of data. Matrix is assumed to have observations in rows
#' @param dist_mthd String specifying method for `dist` function. Default is "euclidian"
#' @param clust_mthd String specifying method for `hclust` function. Default is "ward.D"
#' 
#' @return Value from from `fastcluster::hclust()` called using the arguments
performHierClust <- function(X, dist_mthd="euclidian", clust_mthd="ward.D") {
  fastcluster::hclust(dist(X, method=dist_mthd)^2, method=clust_mthd)
}

#' Compute Pmc for first split in hierarchical clustering for a dataset
#' 
#' Uses the pooled variance estimate to be consistent with the t-test & gao et al.
#' 
#' @param dat N x P data matrix
#' 
#' @return Value of Pmc
computePmcHierClust <- function(dat) {
  clust_n <- length(dat)
  hcl <- performHierClust(dat)
  labels <- cutree(hcl, 2)
  
  ## Get mixture parameters
  mixture_params <- lapply(unique(labels), function(lab) {
    ## Subset the observations
    data_lab <- dat[which(labels == lab)]

    ## Build the parameter object manually
    mcl <- list()
    mcl$d <- 1
    mcl$parameters <- list()
    mcl$G <- 1
    mcl$parameters$mean <- matrix(mean(data_lab, nrow=1))
    mcl$parameters$variance <- list()
    mcl$parameters$variance$sigmasq <- array(var(data_lab), dim=c(1, 1, 1))
    mcl$parameters$pro <- 1
    if (is.null(mcl$parameters$variance$sigmasq)) {
      mcl$parameters$variance$sigmasq <- array(0, dim=c(D, D, mcl$G))
    }
    
    res <- .buildPmcParamsMclust(mcl, T)
    res$prob <- res$prob * length(data_lab) / clust_n
    res$G <- mcl$G
    res
  })
  
  ## Pooled variance estimate
  n1 <- sum(labels == 1)
  n2 <- sum(labels == 2)
  pooled_var <- (n1 - 1) * mixture_params[[1]]$var + (n2 - 1) * mixture_params[[2]]$var
  pooled_var <- pooled_var / (clust_n - 1)
  
  mixture_params[[1]]$var <- pooled_var
  mixture_params[[2]]$var <- pooled_var
  
  computePmc(mixture_params)$integral
}


#' Wrapper around the two-sample t-test for Gao et al.'s procedure
#' 
#' @param dat NxP matrix
#' 
#' @return Output from the `test_hier_clusters_exact() function`
performGaoTest <- function(dat) {
  hcl <- performHierClust(dat)
  test_hier_clusters_exact(dat, link="ward.D", K=2, k1=1, k2=2, hcl=hcl)
}
