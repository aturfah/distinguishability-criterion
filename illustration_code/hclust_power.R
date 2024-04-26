## Results for Section 2.6: Hypothesis testing in hierarchical clustering Power
## Date: 4/20/2024
library(dplyr)
library(tidyr)
library(clusterpval)
library(ggplot2)
library(latex2exp)

rm(list=ls())

source("illustration_code/hyp_test_helpers.R")

## Set the simulation parameters
pmc_cutoff <- 0.094
num_simul <- 500
dist_vec <- seq(1, 7, length.out=13)
sample_size <- 150

set.seed(20240123)

results_raw <- lapply(dist_vec, function(distance) {
  res_d <- sapply(1:num_simul, function(idx) {
    if (idx %% 5 == 0 | idx == 1) cat(paste0("Dist ", distance, "; Simulation #", idx), "\t\r")
    output_filename <- paste0("output/hclust_power/hclust_power_D", distance, "_", idx, ".RData")
    if (file.exists(output_filename)) {
      load(output_filename)
    } else {
      dat <- c(rnorm(sample_size/2), rnorm(sample_size/2, mean=distance))
      dat <- matrix(dat, ncol=1)
      
      pmc <- computePmcHierClust(dat)
      gao <- performGaoTest(dat)
      save(dat, pmc, gao, file=output_filename)
    }

    return(c(pmc=pmc, gao=gao$pval))
  })
  
  t(res_d) %>%
    data.frame(row.names=NULL) %>%
    mutate(pmc_positive=pmc<pmc_cutoff,
           gao_positive=gao<0.05,
           dist=distance)
})


results <- Reduce(rbind, results_raw)


plt <- results %>%
  group_by(dist) %>%
  summarize(pmc_detect=mean(pmc_positive),
            gao_detect=mean(gao_positive)) %>%
  pivot_longer(cols=ends_with("detect"), names_to="Type") %>%
  mutate(Type=factor(Type, levels=c("pmc_detect", "gao_detect"))) %>%
  ggplot(aes(x=dist, y=value, group=Type, linetype=Type)) +
  scale_linetype(
    labels=c(TeX("$P_{MC}$"), "Gao et al."),
    name="Method"
  ) +
  scale_x_continuous(breaks=1:7) +
  geom_point(size=1) + geom_line() +
  ylab("Power") +
  xlab(TeX(r"(|$\mu$ ${}_{1}$ - $\mu$ ${}_{2}$| / $\sigma$)")) +
  theme_bw() +
  theme(text=element_text(size=8),
        legend.margin = margin(0),
        legend.position = "bottom")

plot(plt)

ggsave("plots/hclust_power.png",
       plt,
       height=2, width=6, units="in")
