## 2024/01/29: Take the fpc package functions and make them compatible with the custom functions

library(fpc)
library(class)

predStrength <- function (xdata, Gmin = 2, Gmax = 10, M = 50, clustermethod = kmeansCBI, 
                                     classification = "centroid", centroidname = NULL, cutoff = 0.8, 
                                     nnk = 1, distances = inherits(xdata, "dist"), count = FALSE, 
                                     ...) 
{
  classifnp.mod <- function (data, clustering, method = "centroid", cdist = NULL, 
                             centroids = NULL, nnk = 1) 
  {
    data <- as.matrix(data)
    k <- max(clustering)
    p <- ncol(data)
    n <- nrow(data)
    topredict <- clustering < 0
    if (method == "averagedist") {
      if (is.null(cdist)) {
        cdist <- as.matrix(dist(data))
      }
      else {
        cdist <- as.matrix(cdist)
      }
      prmatrix <- matrix(0, ncol = k, nrow = sum(topredict))
      for (j in 1:k) {
        prmatrix[, j] <- rowMeans(as.matrix(cdist[topredict, 
                                                  clustering == j, drop = FALSE]))
      }
      clpred <- apply(prmatrix, 1, which.min)
      clustering[topredict] <- clpred
    }
    if (method == "centroid") {
      if (is.null(centroids)) {
        centroids <- matrix(0, ncol = p, nrow = k)
        for (j in 1:k) centroids[j, ] <- colMeans(as.matrix(data[clustering == 
                                                                   j, ]))
      }
      k <- nrow(centroids)
      clustering[topredict] <- knn1(centroids, data[topredict, 
      ], 1:k)
    }
    if (method == "qda") {
      qq <- try(qda(data[!topredict, ], grouping = as.factor(clustering[!topredict])), 
                silent = TRUE)
      if (identical(attr(qq, "class"), "try-error")) 
        qq <- lda(data[!topredict, ], grouping = as.factor(clustering[!topredict]))
      clustering[topredict] <- as.integer(predict(qq, data[topredict, 
      ])$class)
    }
    if (method == "lda") {
      qq <- lda(data[!topredict, ], grouping = as.factor(clustering[!topredict]))
      clustering[topredict] <- as.integer(predict(qq, data[topredict, 
      ])$class)
    }
    if (method == "knn") {
      clustering[topredict] <- as.integer(knn(data[!topredict, 
      ], data[topredict, ], as.factor(clustering[!topredict]), 
      k = nnk))
    }
    if (method == "fn") {
      if (is.null(cdist)) {
        cdist <- as.matrix(dist(data))
      }
      else {
        cdist <- as.matrix(cdist)
      }
      fdist <- matrix(0, nrow = sum(topredict), ncol = k)
      for (i in 1:k) {
        fdist[, i] <- apply(as.matrix(cdist[topredict, clustering == 
                                              i]), 1, max)
      }
      bestobs1 <- apply(fdist, 1, which.min)
      clustering[topredict] <- bestobs1
    }
    clustering
  }
  
  
  xdata <- as.matrix(xdata)
  n <- nrow(xdata)
  nf <- c(floor(n/2), n - floor(n/2))
  indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
  corrpred <- list()
  for (k in Gmin:Gmax) {
    if (count) 
      cat(k, " clusters\n")
    corrpred[[k]] <- numeric(0)
    for (l in 1:M) {
      nperm <- sample(n, n)
      if (count) 
        cat(" Run ", l, "\n")
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
      for (i in 1:2) {
        if (distances) 
          clusterings[[i]] <- clustermethod(as.dist(xdata[indvec[[l]][[i]], 
                                                          indvec[[l]][[i]]]), k, ...)
        else clusterings[[i]] <- clustermethod(xdata[indvec[[l]][[i]], 
        ], k, ...)
        jclusterings[[i]] <- rep(-1, n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
        centroids <- NULL
        if (classification == "centroid") {
          if (is.null(centroidname)) {
            if (identical(clustermethod, kmeansCBI)) 
              centroidname <- "centers"
            if (identical(clustermethod, claraCBI)) 
              centroidname <- "medoids"
          }
          if (!is.null(centroidname)) 
            centroids <- clusterings[[i]]$result[centroidname][[1]]
        }
        j <- 3 - i
        if (distances) 
          classifications[[j]] <- classifdist(as.dist(xdata), 
                                              jclusterings[[i]], method = classification, 
                                              centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
        else classifications[[j]] <- classifnp.mod(xdata, 
                                                   jclusterings[[i]], method = classification, 
                                                   centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
      }
      ps <- matrix(0, nrow = 2, ncol = k)
      for (i in 1:2) {
        ctable <- xtable(clusterings[[i]]$partition, 
                         classifications[[i]], k)
        for (kk in 1:k) {
          ps[i, kk] <- sum(ctable[kk, ]^2 - ctable[kk, 
          ])
          cpik <- clusterings[[i]]$partition == kk
          nik <- sum(cpik)
          if (nik > 1) 
            ps[i, kk] <- ps[i, kk]/(nik * (nik - 1))
          else ps[i, kk] <- 1
        }
      }
      corrpred[[k]][l] <- mean(c(min(ps[1, ]), min(ps[2, 
      ])))
    }
  }
  mean.pred <- numeric(0)
  if (Gmin > 1) 
    mean.pred <- c(1)
  if (Gmin > 2) 
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
  for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(corrpred[[k]]))
  optimalk <- max(which(mean.pred > cutoff))
  out <- list(predcorr = corrpred, mean.pred = mean.pred, optimalk = optimalk, 
              cutoff = cutoff, method = clusterings[[1]]$clustermethod, 
              Gmax = Gmax, M = M)
  class(out) <- "predstr"
  out
}
