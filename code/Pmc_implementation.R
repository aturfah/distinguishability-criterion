## Functions that implement the Distinguishability criterion via cubature and Monte Carlo integration
## Date: 4/19/2024

library(cubature) ## Cubature integration package
library(parallel) ## mclapply for MC Integration
library(mvtnorm)  ## Multivariate Gaussian functions (sampling, density)
library(abind)    ## Concatenate arrays


#' Take the output of `Mclust()` and format it for consumption by `compute*Pmc` functions
#' 
#' @param res_mclust Output of `Mclust()` function call
#' @param single_element Whether to combine into a single list element
.buildPmcParamsMclust <- function(res_mclust, single_element=F) {
  params <- res_mclust$parameters
  
  K <- length(params$pro)
  if (res_mclust$d == 1) {
    params$mean <- matrix(params$mean, nrow=1)
    params$variance$sigma <- array(params$variance$sigmasq, dim=c(1, 1, K))
  }
  
  output <- lapply(1:K, function(j) {
    list(
      prob=params$pro[j],
      mean=params$mean[, j, drop=F],
      var=params$variance$sigma[, , j, drop=F]
    )
  })
  
  if (single_element) {
    g <- function(a, b) {
      list(
        mean=cbind(a$mean, b$mean),
        prob=c(a$prob, b$prob),
        var=abind(a$var, b$var, along=3)
      )
    }
    output <- Reduce(g, output)
    
  } else {
    return(output)
  }
}


#' Draw observations from a given Gaussian mixture component.
#' 
#' Generates `num_samples` observations from a Gaussian with parameters described by `distbn_params`
#' 
#' @param distbn_params List containing the mean, variance, and proportions for this mixture component
#' @param num_samples Number of observations to draw from this distribution
#' 
#' @return A matrix with `num_samples` rows
.sampleDistbn <- function(distbn_params, num_samples) {
  mean <- distbn_params$mean
  var <- distbn_params$var
  prob <- distbn_params$prob
  
  if (is.null(nrow(mean))) mean <- matrix(mean, nrow=1, ncol=1)
  
  D <- nrow(mean)
  K <- length(prob)
  
  if (K > 1) {
    ## Split the samples across each sub-distribution
    Nk_vec <- c(rmultinom(1, num_samples, prob))
    
    ## With the counts sample from the subdivision
    output <- lapply(1:K, function(idx) {
      if (Nk_vec[idx] == 0) return(NULL)
      sample_distbn(list(
        mean=mean[, idx, drop=F],
        var=var[, , idx, drop=F],
        prob=1
      ),
      Nk_vec[idx]
      )
    })
    return(Reduce(rbind, output))
  }
  
  rmvnorm(num_samples, mean, matrix(var[, , 1], ncol=D, nrow=D))
}

#' Draw observations from the overall mixture distrubition
#' 
#' Generates `num_samples` observations from the Gaussian mixture distribution with parameters specified in `distbn_params_list`
#' 
#' @param distbn_params_list List of lists containing the mean, variance, and proportions for all mixture components.
#' @param num_samples Total number of samples to draw from the overall mixture distribution.
#' 
#' @return A matrix with `num_samples` rows
.sampleMixture <- function(distbn_params_list, num_samples) {
  prob_vec <- sapply(distbn_params_list, function(x) sum(x$prob))
  Nk_vec <- c(rmultinom(1, num_samples, prob_vec))
  
  distbn_params_list[which(Nk_vec == 0)] <- NULL ## Don't sample components with 0 samples
  Nk_vec <- Nk_vec[which(Nk_vec != 0)]
  
  output <- lapply(1:length(distbn_params_list), function(idx) {
    distbn_params <- distbn_params_list[[idx]]
    .sampleDistbn(distbn_params, Nk_vec[idx])
  })
  
  Reduce(rbind, output)
}

#' PDF of a Gaussian mixture component
#' 
#' For a GMM component with parameters specified by `distbn_params`, return a function to compute its pdf
#' 
#' @param distbn_params List with mean, covariance, and probability measurements
#' 
#' `distbn_params` should have three named components
#' - `mean` should be a `DxK` matrix where each column corresponds to a component mean
#' - `var` should be a `DxDxK` array where each slice corresponds to a covariance matrix
#' - `prob` should be a `K` dimensional vector for the proportions within the parent mixture
.generateDistbnFunc <- function(distbn_params) {
  ## Normalize the probabilities of each component
  M <-  length(distbn_params$prob)
  rel_prop <- distbn_params$prob / sum(distbn_params$prob)  
  
  ## Output function
  function(x) {
    D <- ncol(x)
    res <- sapply(1:M, function(idx) {
      rel_prop[idx] * dmvnorm(x, 
                              distbn_params$mean[, idx, drop=F], 
                              matrix(
                                distbn_params$var[, , idx, drop=F],
                                nrow=D, ncol=D
                              ))
    })
    
    if (is.null(dim(res))) return(sum(res))
    
    return(apply(res, 1, sum))
  }
}

#' Posterior assignment probability ($\pi_j$(x)) function
#' 
#' Get a function to compute the posterior assignment probability for the $j^{th}$ cluster
#' 
#' @param distbn_params_list List containing all parameters for the data mixture distribution. See `generateDistbnFunc` for component structure
#' @param j Index of the cluster for which to generate the function for
#' 
#' @return Function taking argument `x` to compute $\pi_j(x)$ 
.generatePosteriorProbFunc <- function(distbn_params_list, j) {
  function(x) {
    if (is.null(dim(x))) stop("x must be an NxD matrix")
    
    num <- sum(distbn_params_list[[j]]$prob) * .generateDistbnFunc(distbn_params_list[[j]])(x)
    den <- 0
    K <- length(distbn_params_list)
    for (k in 1:K) {
      den <- den + 
        sum(distbn_params_list[[k]]$prob) * .generateDistbnFunc(distbn_params_list[[k]])(x)
    }
    den <- ifelse(den == 0, 1, den) ## Should probably be able to do this in a MC way
    
    return(num / den)
  }
}


#' Monte Carlo $P_{\rm{mc}}$ computation
#'
#' Compute Monte Carlo estimate of Pmc associated with given cluster parameters based on Gaussian cluster assumption.
#' 
#' @param distbn_params_list List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param mc.samples Numeric for number of MC samples to use to approximate the integral. Default 1e5
#' @param obs.per.batch Numeric for the observations to assign to each core. Helps with memory concerns. Default 1e3.
#' @param num.cores Number of cores to use in mclapply call. Default is `floor(detectCores() / 2)`.
#' 
#' @return Value of Pmc from the Monte Carlo integral
computeMonteCarloPmc <- function(distbn_params_list, mc.samples=1e5, obs.per.batch=1e3, num.cores=ceiling(detectCores() / 2)) {
  K <- length(distbn_params_list)
  
  num.batches <- ceiling(mc.samples / obs.per.batch)
  mc.postJ <- mclapply(1:num.batches, function(idx) {
    obs <- .sampleMixture(distbn_params_list, obs.per.batch)
    if (is.null(dim(obs))) obs <- matrix(obs, ncol=1)
    tmp <- sapply(1:K, function(j) {
      postJ <- .generatePosteriorProbFunc(distbn_params_list, j)
      mc.postJ <- postJ(obs)
      sum(mc.postJ * (1 - mc.postJ)) / mc.samples
    })
    sum(tmp)
  }, mc.cores=num.cores)
  
  return(sum(Reduce(c, mc.postJ)))
}


#' Cubature $P_{\rm{mc}}$ computation
#' 
#' Compute Pmc for a given Gaussian mixture distribution based on `cubature` package
#' 
#' @param distbn_params_list List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param integralControl Specifies arguments for integration methods. See details.
#' 
#' @details 
#' TODO: Clear this up
#' For `pcubature` the control variables are `method` for the integration
#' For `pcubature` the control variables are `lowerLimit` and `upperLimit` for the integral (defaults to \pm Inf)
#' For `pcubature` the control variables are `maxEval` which sets the number of integral evaluations (defaults to 1e6)
#' For `pcubature` the control variables are `relTol` which sets the maximum tolerance (defaults to 1e-5)
#' 
#' @return Output from the `cubature::cubintegrate()` function.
computePmc <- function(distbn_params_list, integralControl=list()) {
  
  ## Default Integral Parameters
  intCont <- list(
    method="hcubature",
    lowerLimit=-Inf,
    upperLimit=Inf,
    maxEval=1e6,
    relTol=1e-5,
    nVec=1024L
  )
  intCont[names(integralControl)] <- integralControl
  
  ## Get the parameters
  K <- length(distbn_params_list)
  D <- nrow(distbn_params_list[[1]]$mean)
  
  ## Pmc trivially 0 for a single cluster
  if (K == 1) return(list(
    integral=0,
    error=NULL,
    neval=NULL,
    returnCode=0
  ))
  
  ## Input validation: Make sure component probabilities all scale to 1
  probs <- lapply(distbn_params_list, function(x) x$prob)
  total_prob <- sum(Reduce(c, probs))
  
  if (total_prob != 1) {
    warning("Probabilities in distbn_params_list do not sum to 1. Re-scaling...")
    for (idx in 1:K) {
      distbn_params_list[[idx]]$prob <- distbn_params_list[[idx]]$prob / total_prob
    }
  }
  
  ## Input validation: Make sure dimensions line up for means/covariances
  for (j in 1:K) {
    M <- length(distbn_params_list[[j]]$prob)
    ## Means checking
    stopifnot(D == nrow(distbn_params_list[[j]]$mean))
    stopifnot(M == ncol(distbn_params_list[[j]]$mean))
    
    ## Covariance Matrix checking
    stopifnot(D == dim(distbn_params_list[[j]]$var[1]))
    stopifnot(D == dim(distbn_params_list[[j]]$var[2]))
    stopifnot(M == dim(distbn_params_list[[j]]$var[3]))
  }
  rm(j, M)
  
  ## Functions to perform integral
  pmcIntegralFunc <- function(x, cubatureFunc=T) {
    if (is.null(dim(x))) x <- matrix(x, nrow=1)
    
    if (cubatureFunc == T) x <- t(x)
    
    fX <- sapply(1:K, function(j) {
      comp_prob <- sum(distbn_params_list[[j]]$prob)
      comp_prob * .generateDistbnFunc(distbn_params_list[[j]])(x)
    })
    fX <- apply(fX, 1, sum)
    
    res <- sapply(1:K, function(j) {
      postJFunc <- .generatePosteriorProbFunc(distbn_params_list, j)
      postJ <- postJFunc(x)
      
      res <- postJ * (1 - postJ) 
      if (cubatureFunc == T) {
        res <- res * fX
      }
      return(res)
    })
    if (is.null(dim(res))) return(sum(res))
    
    output <- apply(res, 1, sum)
    if (cubatureFunc == T) return(matrix(output, ncol=nrow(x)))
    
    output
  }
  
  res <- cubintegrate(
    f=pmcIntegralFunc,
    lower=rep(intCont$lowerLimit, D),
    upper=rep(intCont$upperLimit, D),
    nVec=intCont$nVec,
    relTol=intCont$relTol,
    maxEval=intCont$maxEval,
    method=intCont$method
  )
  
  return(res)
}

