## Results for Section 2.6: Hypothesis testing in hierarchical clustering p-values
## Date: 4/20/2024

library(clusterpval)
library(dplyr)
library(tidyr)
library(latex2exp)
library(stringr)
library(ggplot2)
library(ggh4x)

rm(list=ls())


source("code/Pmc_implementation.R")
source("illustration_code/hyp_test_helpers.R")

#' Compute Pmc p-value for first split in hierarchical clustering via bootstrap procedure
#' 
#' @param dat NxP data matrix
#' @param bs_samp Number of bootstrap samples (default is 500)
#' 
#' @return Bootstrap p-value
computePmcBootstrapPvalue <- function(dat, bs_samp=500) {
  orig_pmc <- computePmcHierClust(dat)
  
  ## Resample to get bootstrap p-values
  bs_pmc <- sapply(1:bs_samp, function(idx) {
    if (idx %% 50 == 0 | idx == 1) cat(paste0("Bootstrap Simulation #", idx), "\r")
    n <- length(dat)
    bs_idx <- sample(1:n, n, replace=T)
    computePmcHierClust(dat[bs_idx])
  })
  mean(orig_pmc < bs_pmc)
}

## Set the sample for all data to be evaluated separately
set.seed(20240119)
replicates <- 5000
data <- lapply(1:replicates, function(x) rnorm(150))
output_prefix <- "output/hclust_hyp_test"

if (!dir.exists(output_prefix)) 
  dir.create(output_prefix)

##############################################
##### Pmc Empirical Distribution #############
##############################################

## Compute null distribution and get threshold
set.seed(1)
null_distbn_file <- file.path(output_prefix, "pmc_null_distbn.RData")
if (file.exists(null_distbn_file)) {
  load(null_distbn_file)
} else {
  null_pmc_distbn <- sapply(1:5000, function(idx) {
    if (idx %% 50 == 0 | idx == 1) cat(paste0("Simulation #", idx), "\r")
    
    computePmcHierClust(rnorm(150))
  })
  save(null_pmc_distbn, file=null_distbn_file)
}

pmc_threshold <- unname(quantile(null_pmc_distbn, 0.05))

## Get the Pmc values for the data to be evaluated
emp_distbn_file <- file.path(output_prefix, "pmc_emp_distbn.RData")
if (file.exists(emp_distbn_file)) {
  load(emp_distbn_file)
} else {
  emp_pmc <- sapply(1:length(data), function(idx) {
    if (idx %% 50 == 0 | idx == 1) cat(paste0("Simulation #", idx), "\r")
    dat <- data[[idx]]
    computePmcHierClust(dat)
  })
  save(emp_pmc, file=emp_distbn_file)
}
rm(null_distbn_file, emp_distbn_file)

## Compute the p-values
empirical_pvalue <- sapply(emp_pmc, function(x) mean(x < null_pmc_distbn))

##############################################
##### Pmc Bootstrap Distribution #############
##############################################

bootstrap_pvalue <- sapply(1:length(data), function(idx) {
  ## These simulations take a hot minute so should write them to disk
  output_file <- file.path(output_prefix, paste0("bootstrap_", idx, ".Rdata"))
  if (file.exists(output_file)) {
    load(output_file)
  } else {
    output <- computePmcBootstrapPvalue(data[[idx]])
    save(output, file=output_file)
  }
  if (idx == 1 | idx %% 5 == 0)  cat("Completed:", idx, "    \n\n")
  
  return(output)
})

###########################################
########## Gao et al. p-value #############
###########################################

gao_pvalue <- sapply(1:length(data), function(idx) {
  if (idx %% 10 == 0 | idx == 1) cat(paste0("Gao Simulation #", idx), "\r")
  
  output_file <- file.path(output_prefix, paste0("gao_", idx, ".Rdata"))
  if (file.exists(output_file)) {
    load(output_file)
  } else {
    dat <- as.matrix(data[[idx]], ncol=1)
    output <- performGaoTest(dat)
    save(output, file=output_file)
  }
  return(output$pval)
})


###########################################
########## 2 Sample T p-value #############
###########################################

ttest_pvalue <- sapply(1:length(data), function (idx) {
  if (idx %% 100 == 0 | idx == 1) cat(paste0("Simulation #", idx), "\r")
  
  dat <- data[[idx]]
  hcl <- performHierClust(dat)
  labs <- cutree(hcl, 2)
  
  X1 <- dat[which(labs==1)]
  X2 <- dat[which(labs==2)]
  
  res <- t.test(X1, X2, alternative="two.sided", var.equal=T)
  res$p.value
})


###############################################
####### Generate results and figures ##########
###############################################

print(paste("Empirical Pmc 0.05 threshold:", round(pmc_threshold, 3)))
print(paste("Empirical Pmc Type I Error Rate:", mean(emp_pmc < round(pmc_threshold, 3))))
print(paste("Bootstrap Pmc Type I Error Rate", mean(bootstrap_pvalue < 0.05)))

## Combine these to get plot
results <- data.frame(
    Pmc=empirical_pvalue,
    PmcBS=bootstrap_pvalue,
    Gao=gao_pvalue,
    ttest=ttest_pvalue
  ) %>%
  pivot_longer(cols=c("Pmc", "PmcBS", "Gao", "ttest"), names_to="mthd") %>%
  mutate(mthd=ifelse(mthd=="ttest", "2SampleT", mthd),
         mthd=factor(mthd, levels=c("Pmc", "PmcBS", "Gao", "2SampleT")))

facet_labels <- c("Pmc"=TeX("$P_{MC}$", output="character"),
                  "PmcBS"=TeX("$P_{MC}$ Bootstrap", output="character"),
                  "Gao"="Gao~et~al.",
                  "2SampleT"=TeX("Two~Sample~t-Test"))

levels(results$mthd) <- facet_labels

plt <- results %>%
  filter(str_count(mthd, "Bootstrap") == 0) %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=25,
                 color="black",
                 fill="white",
                 boundary=c(0),
                 closed="right") +
  facet_wrap(. ~ mthd, scales="free",
             labeller=label_parsed
  ) +
  facetted_pos_scales(
    x = list(
      mthd != "P[MC]" ~ scale_x_continuous(limits=c(0, 1))
    ),
    y = list(
      mthd == "Two~Sample~t-Test" ~ scale_y_continuous(limits=c(0, 5100)),
      mthd == "P[MC]*' Bootstrap'" ~ scale_y_continuous(limits=c(0, 600)),
      mthd == "P[MC]" ~ scale_y_continuous(limits=c(0, 600)),
      mthd == "Gao~et~al." ~ scale_y_continuous(limits=c(0, 600))
    )
  ) +
  xlab("") + ylab("Count") +
  theme_bw() + theme(text=element_text(size=8))

plot(plt)

ggsave("plots/hclust_hyp_test.png", 
       plot=plt,
       width=6, height=1.5, units="in")


