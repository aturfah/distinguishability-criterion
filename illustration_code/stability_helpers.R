## Functions to perform clustering and compute cluster stability measures
## Date: 04/26/2024


library(ClusterR)
library(cluster)
library(mclust)
library(pdfCluster)
library(class)
library(parallel)

source("code/Pmc_implementation.R")
source("illustration_code/prediction_strength_function.R")
source("illustration_code/hyp_test_helpers.R")

#' Function to perform k-means clustering
#' 
#' @param dat NxP matrix of data
#' @param k Number of clusters to fit
#' 
#' @return Output from `ClusterR::KMeans_rcpp()`
kmeansFunc <- function(dat, k) {
  if (is.null(dim(dat))) dat <- matrix(dat, ncol=1)
  res <- KMeans_rcpp(dat, k, num_init=5)
  res$cluster <- res$clusters
  
  res$cluster <- res$clusters
  res$partition <- res$clusters
  res$nc <- k
  res$clusterlist <- lapply(1:k, function(lab) res$clusters == lab) ## List of points belonging to each cluster
  res$clusterMethod <- "kmeansCustom"
  
  res$result <- list()
  res$result$centroid <- res$centroids
  
  res
}


#' Function to perform hierarchical clustering
#' @param dat NxP data matrix
#' @param k Number of clusters to split
#' 
#' @return Object in the same format as `ClusterR::KMeans_rcpp()`
hclFunc <- function(dat, k) {
  hcl <- performHierClust(dat)
  # partitions <- min_size_cluster_partitions(hcl, min_size=1, max_clusters=max(K_VEC))
  cluster_labels <- cutree(hcl, k)
  centroids <- sapply(1:k, function(lab) {
    dat_clust <- dat[which(cluster_labels == lab), , drop=F]
    colMeans(dat_clust)
  })
  list(cluster=cluster_labels,
       clusters=cluster_labels,
       partition=cluster_labels,
       centroid=t(centroids),
       centroids=t(centroids))
}


#' Function to compute cluster stability based on Lange et al.'s method using ARI as the similarity measure
#' 
#' @param dat NxP matrix of data
#' @param k Number of clusters
#' @param clustFunc Function to perform clustering; arguments are dat and k
#' @param B Number of times to resample to compute stability
#' @param num_cores Number of cores for `mclapply`
#' 
#' @return Numeric value for stability measure
stabilityARI <- function(dat, k, clustFunc, B=100, num_cores=1) {
  if (is.null(dim(dat))) data <- matrix(dat, ncol=1)
  N <- nrow(dat)
  
  stability <- mclapply(1:B, function(idx) {
    set.seed(idx)
    
    splitY <- sample(1:N, N/2)
    splitZ <- sample(1:N, N/2)
    
    datY <- dat[splitY, , drop=F]
    datZ <- dat[splitZ, , drop=F]

    ## Get clusterings of Y and Z
    clusterY <- clustFunc(datY, k)
    labelsY <- clusterY$cluster
    clusterZ <- clustFunc(datZ, k)
    labelsZ <- clusterZ$cluster
    
    ## Predict clustering of Z based on nearest point in Y
    labelsZPred <- class::knn(train=datY, cl=labelsY, test=datZ)
    
    ## Also compare ARI as cluster agreement
    ari <- pdfCluster::adj.rand.index(labelsZPred, labelsZ)
    
    return(ari)
  }, mc.cores=num_cores)
  
  mean(unlist(stability))
}


#' Function to compute Stability Measures + Pmc for a given clustering method
#' 
#' @param dat NxP matrix of data
#' @param k_values Range of number of clusters to search
#' @param clustFunc Function to perform clustering; arguments are dat and k
#' 
#' @return Data frame for the stability measures table
clusterResults <- function(dat, k_values, clustFunc=kmeansFunc) {
  ## Compute the cluster-specific terms
  results_list <- lapply(k_values, function(k) {
    res <- clustFunc(dat, k)
    
    data_labels <- res$cluster
    N <- length(data_labels)
    
    ## Get the cluster parameters
    clust_params <- lapply(1:k, function(lab) {
      clust_idx <- which(data_labels == lab)
      
      clust_dat <- dat[clust_idx, , drop=F]
      
      if (nrow(clust_dat) == 1) {
        D <- ncol(clust_dat)
        params <- list(prob=1,
                       mean=t(clust_dat),
                       var=array(0, dim=c(D, D, 1)))
      } else {
        mcl <- Mclust(clust_dat, G=1, verbose=F)
        params <- .buildPmcParamsMclust(mcl, T)        
      }

      params$prob <- mean(data_labels == lab)
      params
    })

    ## Output
    c(k=k,
      Pmc=computePmc(clust_params)$integral,
      Silhouette=silhouette_of_clusters(dat, data_labels)$silhouette_global_average,
      Stability=stabilityARI(dat, k, clustFunc))
  })
  
  ## Build data frame
  results <- data.frame(Reduce(rbind, results_list),
                        row.names=NULL)
  
  ## Add Gap Statistic + Prediction Strength
  gap_results <- clusGap(dat, clustFunc, K.max=max(k_values),
                         B=500, d.power=2,
                         verbose=F)
  pred_str <- predStrength(dat, Gmax=max(k_values),
                           M=100,
                           clustermethod=kmeansFunc,
                           centroidname="centroid")
  
  results$Gap <- gap_results$Tab[, "gap"]
  results$Pred.Strength <- pred_str$mean.pred
  
  return(results)
}

