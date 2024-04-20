## Results for Section 4.5: Monte Carlo vs Cubature Pmc evaluation

library(dplyr)

rm(list=ls())

source("code/Pmc_implementation.R")

#' Generate parameters for the simulation
#' 
#' @param distance Euclidean distance between cluster centers
#' @param dim Dimension for the distributions
#' 
#' @return List of MVN parameters
generateDistbnParams <- function(distance, dim) {
  simulParams <- list(
    list(
      prob=1/3,
      mean=matrix(sqrt(distance^2/dim), nrow=dim),
      var=array(diag(nrow=dim), dim=c(dim, dim, 1))
    ),
    list(
      prob=1/3,
      mean=matrix(0, nrow=dim),
      var=array(diag(nrow=dim), dim=c(dim, dim, 1))
    ),
    list(
      prob=1/3,
      mean=matrix(-sqrt(distance^2/dim), nrow=dim),
      var=array(diag(nrow=dim), dim=c(dim, dim, 1))
    )
  )
}

replicates <- 50
results <- lapply(1:4, function(dim) {
  print(paste("Dim:", dim))
  
  output_file <- paste0("output/mc_comparison/mc_simul_dim", dim, ".RData") 
  
  if (file.exists(output_file)) {
    cat("\tLoaded!\n")
    load(output_file)
  } else {
    simul_par <- generateDistbnParams(3, dim)
    
    timeMC <- system.time({
      resMC <- sapply(1:replicates, function(x) computeMonteCarloPmc(simul_par, num.cores=8))
    })["elapsed"]
    cat("\tMC:", unname(timeMC)/replicates, "\n")
    
    time <- system.time({
      res <- lapply(1:replicates, function(x) computePmc(simul_par))
    })["elapsed"]
    cat("\tNon-MC:", unname(time)/replicates, "\n")
    
    output <- list(dim=dim,
                   replicates=replicates,
                   res_reg=res,
                   time_reg=unname(time),
                   res_mc=resMC,
                   time_mc=unname(timeMC))
    save(output, file=output_file)
  }
  return(output)
})


results_table <- lapply(results, function(datum) {
  res_reg <- datum$res_reg
  
  pmis_reg <- mean(sapply(res_reg, function(x) x$integral))
  err_reg <- mean(sapply(res_reg, function(x) x$error))
  sd_reg <- sd(sapply(res_reg, function(x) x$integral))
  
  
  res_mc <- datum$res_mc
  pmis_mc <- mean(res_mc)
  sd_mc <- sd(res_mc)
  data.frame(
    D=datum$dim,
    Repl=datum$replicates,
    pmis_reg=pmis_reg,
    err_reg=err_reg,
    sd_reg=sd_reg,
    time_reg=datum$time_reg / datum$replicates,
    pmis_mc=pmis_mc,
    sd_mc=sd_mc,
    time_mc=datum$time_mc / datum$replicates
  )
})

Reduce(rbind, results_table) %>%
  select(D, pmis_reg, time_reg, pmis_mc, sd_mc, time_mc) %>%
  xtable::xtable(digits=5) %>%
  xtable::print.xtable(include.rownames=F)




