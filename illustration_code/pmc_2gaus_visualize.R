## Generate Supplementary Figure 1: Two Gaussians at different levels of Pmc
## NOTE: Should be run from base directory
## Date: 4/19/2024

library(tidyr)
library(dplyr)
library(ggplot2)
library(latex2exp)

rm(list=ls())

source("Pmc_implementation.R")

#' Generate the distribution parameter object for this simulation
#' 
#' @param mu2 Mean for Gaussian distribution not centered at 0
#' 
#' @return List containing parameters for consumption by `computePmc()`
generateDistbnParams <- function(mu2) {
  list(
    list(
      prob=0.5,
      mean=matrix(0, nrow=1, ncol=1),
      var=array(1, dim=c(1, 1, 1))
    ),
    list(
      prob=0.5,
      mean=matrix(mu2, nrow=1, ncol=1),
      var=array(1, dim=c(1, 1, 1))
    )
  )
}

## Cutoff values for which we would like to plot
cutoffs <- c(0.49, 0.4, 0.3, 0.1, 0.05, 0.01)
means <- numeric(length(cutoffs))
realized_pmc <- numeric(length(cutoffs))
plots <- vector("list", length(cutoffs))

## Slowly pull means apart until certain Pmc is reached; store these means
mu2 <- 0
pmc <- computePmc(generateDistbnParams(mu2))$integral
step = 1e-3
for (idx in 1:length(cutoffs)) {
  while(pmc > cutoffs[idx]) {
    mu2 <- mu2 + step
    
    pmc <- computePmc(generateDistbnParams(mu2))$integral
  }
  means[idx] <- mu2
  realized_pmc[idx] <- pmc
  print(paste("Mean for Pmc =", round(pmc, 4), "@", round(mu2, 4)))
  
  ## Visualize this distribution
  distbn_params <- generateDistbnParams(mu2)
  xvals <- seq(-4, 11, length.out=1000)
  dens1 <- dnorm(xvals, mean=distbn_params[[1]]$mean, sd=sqrt(distbn_params[[1]]$var))
  dens2 <- dnorm(xvals, mean=distbn_params[[2]]$mean, sd=sqrt(distbn_params[[2]]$var))
  
  plt <- data.frame(x=xvals, clust1=dens1, clust2=dens2) %>%
    pivot_longer(cols=starts_with("clust")) %>%
    ggplot(aes(x=x, y=value, group=name, linetype=name)) +
    geom_line() +
    xlab("") + ylab("") +
    ggtitle(TeX(paste0("$P_{MC} = ", cutoffs[idx], "$"))) +
    theme_bw() + theme(text=element_text(size=8), legend.position="none")
  
  plots[[idx]] <- plt
}
## Keep the environment clean
rm(idx, mu2, pmc, plt, dens1, dens2, xvals)

## Make sure these are very close
data.frame(
    truth=cutoffs,
    value=realized_pmc,
    diff=cutoffs-realized_pmc
  ) %>%
  t()

## Visualize the distributions
plot_arrange <- do.call(gridExtra::arrangeGrob, c(list(ncol=3, nrow=2), plots))
plot(plot_arrange)

ggsave(filename="plots/pmc_2gaus_visualize.png",
       plot=plot_arrange,
       width=6, height=3, units="in")