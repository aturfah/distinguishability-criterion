## Results for Section 2.5: K-means and constrained optimization procedure
## Date: 04/26/2024

library(mvtnorm)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(latex2exp)

rm(list=ls())

source("code/Pmc_implementation.R")
source("illustration_code/stability_helpers.R")

set.seed(20240124)
dat <- rbind(
  mvtnorm::rmvnorm(150, c(0, 0), diag(nrow=2)),
  mvtnorm::rmvnorm(150, c(-4, 4), diag(nrow=2)),
  mvtnorm::rmvnorm(150, c(1.75, 1.75), diag(nrow=2))
)
true_labels <- c(rep(1, 150), rep(2, 150), rep(3, 150))


## Produce the table with the k-means stability measures
results_kmeans <- clusterResults(dat, 1:7, kmeansFunc)

## Produce the visualization
plot_df <- data.frame(dat, label=true_labels,
                      pred=kmeansFunc(dat, 3)$clusters)

plt_data <- plot_df %>%
  ggplot(aes(x=X1, y=X2, fill=as.factor(label), 
             shape=as.factor(pred))) +
  geom_point(size=0.8) + xlab("") + ylab("") +
  scale_shape_manual(values=c(`3`=21, `2`=22, `1`=24)) +
  theme_bw() + theme(text=element_text(size=8),
                     legend.position = "none")

plt_pmc <- results_kmeans %>%
  cbind(ZERO=1e-6) %>%
  ggplot(aes(x=k, y=ZERO)) +
  geom_segment(aes(x=k, xend=k, y=ZERO, yend=Pmc), lineend="round", linewidth=1.5) +
  ylab(TeX("$P_{MC}$")) +
  scale_x_continuous(breaks=results_kmeans$k) + 
  theme_bw() + theme(text=element_text(size=8), legend.position="none",
                     panel.grid.minor=element_blank())

plt_gap <- results_kmeans %>%
  ggplot(aes(x=k, y=Gap)) +
  scale_x_continuous(breaks=results_kmeans$k) + 
  ylab("Gap Statistic") +
  geom_point(size=0.75) + geom_line() +
  theme_bw() + theme(text=element_text(size=8),
                     panel.grid.minor=element_blank())

plt_gap <- results_kmeans %>%
  cbind(ZERO=1e-6) %>%
  ggplot(aes(x=k, y=ZERO)) +
  geom_segment(aes(x=k, xend=k, y=ZERO, yend=Gap), lineend="round", linewidth=1.5) +
  ylab(TeX("Gap Statistic")) +
  scale_x_continuous(breaks=results_kmeans$k) + 
  theme_bw() + theme(text=element_text(size=8), legend.position="none",
                     panel.grid.minor=element_blank())

plt_silhouette <- results_kmeans %>%
  cbind(ZERO=1e-6) %>%
  mutate(Silhouette=ifelse(k==1, 0, Silhouette)) %>%
  ggplot(aes(x=k, y=ZERO)) +
  geom_segment(aes(x=k, xend=k, y=ZERO, yend=Silhouette), lineend="round", linewidth=1.5) +
  ylab(TeX("Silhouette")) +
  scale_x_continuous(breaks=results_kmeans$k) + 
  theme_bw() + theme(text=element_text(size=8), legend.position="none",
                     panel.grid.minor=element_blank())

combined_plot_nodat <- arrangeGrob(
  plt_gap, plt_pmc, plt_silhouette,
  nrow=1,
  widths=list( unit(1.75, "in"), unit(1.75, "in"), unit(1.75, "in"))
)

combined_plot <- arrangeGrob(
  plt_data, combined_plot_nodat, nrow=2,
  heights=list(unit(1.75, "in"), unit(1.25, "in"))
)

plot(combined_plot)
ggsave("plots/kmeans_example.png",
       combined_plot,
       height=3, width=6, units="in")

## Print the table
results_kmeans %>%
  relocate(Gap, .after=k) %>%
  knitr::kable(digits=4) %>%
  print()
