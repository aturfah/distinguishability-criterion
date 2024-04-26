## Results for Section 2.7.1: Analysis of Palmer Penguins data
## Date: 04/26/2024

rm(list=ls())

library(palmerpenguins)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(latex2exp)

rm(list=ls())

source("code/Pmc_implementation.R")
source("illustration_code/stability_helpers.R")

dat_raw <- penguins[complete.cases(penguins[, c("species", "bill_length_mm", "flipper_length_mm", "sex")]), ] %>%
  filter(sex=="female")

dat <- apply(dat_raw[, c(3, 5)], 2, function(vec) {
  (vec - mean(vec)) / sd(vec)
})


## k-means and hierarchical clustering analysis
results_kmeans <- clusterResults(dat, 1:8, kmeansFunc)
results_hclust <- clusterResults(dat, 1:8, hclFunc)

results_kmeans %>%
  relocate(Gap, .after=k) %>%
  knitr::kable(digits=4) %>%
  print()

results_hclust %>%
  relocate(Gap, .after=k) %>%
  knitr::kable(digits=4) %>%
  print()
