## Results for Section 2.7.2: Inferring Population Structure from HGDP data
## Date: 4/20/2024

library(stringr)

rm(list=ls())
source("code/PHM_algorithm.R")
source("code/visualizations.R")


data_raw <- read.table("../hierarchical-clustering-distinguishability-criterion/real_data/hgdp.txt") %>%
  data.frame()

allele_dosage_matrix <- apply(data_raw[, 4:ncol(data_raw)], 2, function(vec) {
  ref_allele <- substr(vec[1], 1, 1) 
  
  str_count(vec, ref_allele)
})


## PCA
pc_res <- prcomp(allele_dosage_matrix)
pc_matrix <- pc_res$x


## Screeplot
max_pcs <- 15
screeplot <- data.frame(
    x=1:max_pcs,
    y=pc_res$sdev[1:max_pcs]
  ) %>%
  ggplot(aes(x=x, y=y)) +
  geom_point() + geom_line() +
  ylab("Standard Deviation") + xlab("") +
  scale_x_continuous(breaks=1:max_pcs) +
  theme_bw() + theme(panel.grid.minor=element_blank())

plot(screeplot)

ggsave("plots/hgdp_screeplot.png",
       plot=screeplot,
       height=3, width=6, units="in")

## PHM Algorithm
pc_matrix <- pc_matrix[, 1:5]

set.seed(20240420)
gmm_hgdp <- Mclust(pc_matrix)

gmm_hgdp$G

phm_res <- PHM(gmm_hgdp, mc.samples=1e6, num.cores=7)

plotPHMDendrogram(phm_res)
