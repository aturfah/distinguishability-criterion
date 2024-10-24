## Functions that implement the PHM Algorithm
## Date: 4/19/2024

library(mclust)

source("code/Pmc_implementation.R")

#' Monte Carlo $\Delta P_{\rm {mc}}$ Matrix computation
#' 
#' Compute the merging Pmc reduction matrix estimated via Monte Carlo integration
#' 
#' See `computeMonteCarloPmc` for a breakdown of the parameters
computeMonteCarloDeltaPmcMatrix <- function(distbn_params_list, mc.samples=1e5, obs.per.batch=1e3, num.cores=ceiling(detectCores() / 2)) {
  K <- length(distbn_params_list)
  
  output <- matrix(0, K, K)
  num.batches <- ceiling(mc.samples / obs.per.batch)
  # mc.obs <- .sampleMixture(distbn_params_list = , mc.samples)
  cat("MC Samples\n")
  mc.obs <- Reduce(rbind, 
                  mclapply(1:num.batches, function(idx) {
                    lb <- 1 + obs.per.batch * (idx - 1)
                    ub <- min(obs.per.batch * idx, mc.samples)
                    out <- .sampleMixture(distbn_params_list, (ub-lb))
                    if (is.null(dim(out))) out <- matrix(out, ncol=1)

                    return(out)
                  }, mc.cores=num.cores))

  cat("Computing Posterior Matrix\n")
  distbn.mat <- Reduce(cbind, 
            mclapply(1:K, function(j) {
              densJ <- .generateDistbnFunc(distbn_params_list[[j]])
              sum(distbn_params_list[[j]]$prob) * densJ(mc.obs)
            }, mc.cores=num.cores))
  post.mat <- t(apply(distbn.mat, 1, function(x) x / sum(x)))

  cat("Evaluating Integral\n")
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      post.i <- post.mat[, i]
      post.j <- post.mat[, j]

      output[i, j] <- sum(post.i * post.j) / mc.samples
      output[j, i] <- output[i, j]
    }
  }

  output
}

#'  $\Delta P_{\rm {mc}}$ Matrix computation
#' 
#' Compute the merging Pmc reduction matrix calculated using `cubature` package
#' 
#' See `computePmc` for a breakdown of the parameters
computeDeltaPmcMatrix <- function(distbn_params_list, integralControl=list()) {
  K <- length(distbn_params_list)
  D <- nrow(distbn_params_list[[1]]$mean)
  output <- matrix(0, nrow=K, ncol=K)
  
  intCont <- list(
    method="cuhre",
    lowerLimit=-Inf,
    upperLimit=Inf,
    maxEval=1e6,
    relTol=1e-6,
    nVec=1024L
  )
  intCont[names(integralControl)] <- integralControl
  
  
  generatePmcReductionFunc <- function(i, j) {
    post_i <- .generatePosteriorProbFunc(distbn_params_list, i)
    post_j <- .generatePosteriorProbFunc(distbn_params_list, j)
    
    function (x) {
      if (is.null(dim(x))) x <- matrix(x, nrow=1)
      x <- t(x) ## Cubature passes in transposed to what we expect
      
      fX <- sapply(1:K, function(k) {
        comp_prob <- sum(distbn_params_list[[k]]$prob)
        comp_prob * .generateDistbnFunc(distbn_params_list[[k]])(x)
      })
      fX <- apply(fX, 1, sum)
      
      output <- post_i(x) * post_j(x) * fX
      return(matrix(output, ncol=nrow(x)))
    }
  }
  
  for (i in 1:(K-1)) {
    for (j in ((i+1):K)) {
      res <- cubintegrate(
        f=generatePmcReductionFunc(i, j),
        lower=rep(intCont$lowerLimit, D),
        upper=rep(intCont$upperLimit, D),
        nVec=intCont$nVec,
        method=intCont$method,
        maxEval=intCont$maxEval,
        relTol=intCont$relTol
      )
      
      output[i, j] <- res$integral
      output[j, i] <- output[i, j]
    }
  }
  
  output
}

#' Data posterior matrix computation
#' 
#' Given a mixture-of-mixtures distribution, compute the posterior cluster assignment probability
#' 
#' @param distbn_params_list List of lists of mixture component parameters
#' @param data NxD matrix of observations.
.computePosteriorProbMatrix <- function(distbn_params_list, data) {
  if (is.null(dim(data))) data <- matrix(data, ncol=1)
  K <- length(distbn_params_list)
  distbn.mat <- Reduce(cbind, 
              lapply(1:K, function(j) {
                densJ <- .generateDistbnFunc(distbn_params_list[[j]])
                sum(distbn_params_list[[j]]$prob) * densJ(data)
              }))
  t(apply(distbn.mat, 1, function(x) x / sum(x)))
}


#' PHM Algorithm
#' 
#' Implements Baudry et al. 2010 merging procedure based on Pmc
#' 
#' @param gmm_res Output of `Mclust` for estimating a GMM
#' @param data NxD matrix of observations. If provided will compute labels based on posterior assignment probability.
#' @param mc_est Whether to use MC-Approximation for Pmc calculations
#' @param ... Additional parameters for the `computePmc*Matrix` function
#' 
#' TODO: Fill in Returns + More detailed description
#' 
#' @return FILL ME IN
PHM <- function(res_mclust, data=NULL, mc_est=T, paramsList=NULL, ...) {
  ## Construct components necessary for Pmc estimation
  if (is.null(paramsList)) {
    paramsList <- .buildPmcParamsMclust(res_mclust)
    K <- res_mclust$G
  } else {
    if (!is.null(res_mclust)) warn("Ignoring provided res_mclust for parameters")
    K <- length(paramsList)
  }

  ## Construct the posterior data labels
  cat("Constructing Posterior Data Matrix\n")
  if (!is.null(data) | is.null(res_mclust)) {
    posterior_mat <- .computePosteriorProbMatrix(paramsList, data)
    data_labels <- apply(posterior_mat, 1, which.max)
    N <- nrow(data)
  } else {
    data_labels <- res_mclust$classification
    N <- res_mclust$n
  }
  
  ## Compute Pmc matrix
  cat("Constructing Pmc Matrix\n")
  pmc_red_mat <- if (mc_est) {
    computeMonteCarloDeltaPmcMatrix(paramsList, ...)
  } else {
    computeDeltaPmcMatrix(paramsList, ...)
  }
  ## Zero out diagonal elements
  diag(pmc_red_mat) <- 0
  
  ## Parameters/data to store
  output <- lapply(1:K, function(k) list(
    clusters=k,
    posterior_matrix=if(is.null(data)) {NULL} else {matrix(1, nrow=N)},
    labels=rep(1, N),
    pmc_change=NA,
    merge_components=c(-1, -1),
    pmc=0))
  
  pmc <- sum(pmc_red_mat)
  TEMP <- pmc_red_mat
  
  for (idx in K:2) {
    ## Store the results
    output[[idx]]$labels <- data_labels
    output[[idx]]$pmc <- sum(TEMP)
    if (idx < K) {
      output[[idx]]$pmc_change <- output[[idx+1]]$pmc - output[[idx]]$pmc
    }
    output[[idx]]$pmc_matrix <- TEMP
    
    
    if (!is.null(data)) {
      output[[idx]]$posterior_matrix <- posterior_mat
    }
    
    ## Identify the components to merge
    maxval <- max(TEMP, na.rm=T)
    cand_rows <- which(TEMP == maxval, arr.ind=T)
    cand_rows <- cand_rows[which(cand_rows[, 1] < cand_rows[, 2]), , drop=F]
    i <- cand_rows[1, 1]
    j <- cand_rows[1, 2]
    output[[idx]]$merge_components <- c(i, j) ## store this
    
    ## Calculate the merged row values
    row_i <- TEMP[i, ]
    row_j <- TEMP[j, ]
    
    new_row <- row_i + row_j
    new_row[i] <- 0
    new_row <- new_row[-j]
    
    ## Shrink matrix my removing j row/column
    TEMP <- TEMP[-j, , drop=F]
    TEMP <- TEMP[, -j, drop=F]
    
    ## Replace i terms with i+j
    TEMP[i, ] <- new_row
    TEMP[, i] <- new_row
    
    ## Update the cluster labels; make sure always 1:K
    if (is.null(data)) {
      data_labels[which(data_labels == j)] <- i
      cluster_ids <- sort(unique(data_labels))
      data_labels <- sapply(data_labels, function(lab) {
        which(cluster_ids == lab)
      })      
    } else {
      posterior_mat[, i] <- posterior_mat[, i] + posterior_mat[, j]
      posterior_mat <- posterior_mat[, -j, drop=F]
      data_labels <- apply(posterior_mat, 1, which.max)
    }
  }
  
  output[[1]]$pmc_change <- output[[2]]$pmc
  
  output  
}


#' Find the solution for PHM at a specified threshold
#' 
#' @param phm_output Output of `PHM`
#' @param threshold Pmc threshold, default is 0.01
#' 
#' TODO: Fill in Returns + More detailed description
#' 
#' @return Cluster results based on GMM satisfying specified threshold
thresholdPHM <- function(phm_output, threshold=0.01) {
  kappa <- length(phm_output)
  for (k in kappa:1) {
    if (phm_output[[k]]$pmc < threshold) {
      break 
    }
  }
  
  return(phm_output[[k]])
}
