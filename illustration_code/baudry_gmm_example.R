## Results for Section 2.4: Simulated data from Baudry et al.
## Date: 4/20/2024

library(mclust)
library(rhdf5)
library(gridExtra)

rm(list=ls())

source("code/PHM_algorithm.R")
source("code/visualizations.R")

data <- h5read("../hierarchical-clustering-distinguishability-criterion/BaudryData/4.1.mat", "/exp/data")[[1]]
colors <- c(`1`="coral4", `2`="royalblue", `3`="coral", `4`="olivedrab", `5`="darkgray", `6`="plum")

## Run the PHM algorithm to get the merges
model <- Mclust(data)
cluster_labels <- model$classification
phm_output <- PHM(model, data=data, mc_est=F)

## Visualizations
dendro <- plotPHMDendrogram(phm_output, colors=colors)
htmp <- plotPmcMatrix(phm_output, colors=colors, include_pmc_title = F)
scatterplot <- data %>%
  data.frame(Cluster=factor(cluster_labels)) %>%
  ggplot(aes(x=X1, y=X2, fill=Cluster)) +
  geom_point(size=1, shape=21) +
  scale_fill_manual(values=colors) + 
  xlab("") + ylab("") +
  theme_bw() + theme(text=element_text(size=8),
                     legend.position = "none") 

combined <- arrangeGrob(scatterplot, htmp, 
                        nrow=1, ncol=2, 
                        widths=list( unit(3, "in"), unit(3, "in") ))

plot(combined)

ggsave("plots/baudry_simulation_visualize.png",
       combined,
       units="in", height=2.5, width=6)
