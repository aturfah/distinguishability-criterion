## Results for Section 2.7.2: Inferring Population Structure from HGDP data
## Date: 4/20/2024

library(stringr)
library(gridExtra)

rm(list=ls())
source("code/PHM_algorithm.R")
source("code/visualizations.R")


data_raw <- read.table("../hierarchical-clustering-distinguishability-criterion/real_data/hgdp.txt") %>%
  data.frame()

allele_dosage_matrix <- apply(data_raw[, 4:ncol(data_raw)], 2, function(vec) {
  ref_allele <- substr(vec[1], 1, 1) 
  
  str_count(vec, ref_allele)
})
geographic_labels <- factor(data_raw[, 3],
                            levels=c("EUROPE",
                                     "CENTRAL_SOUTH_ASIA",
                                     "MIDDLE_EAST", 
                                     "OCEANIA",
                                     "EAST_ASIA",
                                     "AMERICA",
                                     "AFRICA"), ordered=T)
geographic_levels <- levels(geographic_labels)
geographic_levels <- case_when(
  geographic_levels == "CENTRAL_SOUTH_ASIA" ~ "C/S Asia",
  geographic_levels == "MIDDLE_EAST" ~ "Middle East",
  geographic_levels == "EUROPE" ~ "Europe",
  geographic_levels == "AFRICA" ~ "Africa",
  geographic_levels == "OCEANIA" ~ "Oceania",
  geographic_levels == "EAST_ASIA" ~ "East Asia",
  geographic_levels == "AMERICA" ~ "America",
  TRUE ~ geographic_levels
)
levels(geographic_labels) <- geographic_levels

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

phm_res <- PHM(gmm_hgdp, data=pc_matrix, mc.samples=1e6, num.cores=7)


colors <- c("posterior.1"="coral",
            "posterior.2"="olivedrab",
            "posterior.3"="plum",
            "posterior.4"="snow4", 
            "posterior.5"="gold",
            "posterior.6"="coral4",
            "posterior.7"="royalblue")


plt_dendro <- plotPHMDendrogram(phm_res, colors=colors)
plt_distruct <- plotPHMDistruct(phm_res, labels=geographic_labels, colors=colors,
                                include_title=T)

plt_dendro
plt_distruct

plt_joined <- arrangeGrob(plt_dendro, 
                          plt_distruct, 
                          nrow=2, ncol=1,
                          heights=list( unit(3, "in"), unit(1, "in")  ))
plot(plt_joined)
# ggsave("meeting_notes/2024_0214/hgdp_results/dendogram_distruct.png",
#        dend_disrupt_joined,
#        width=6, height=4, units="in")
ggsave("plots/hgdp_phm.png",
       plot=plt_joined,
       width=6, height=4,
       units="in")


