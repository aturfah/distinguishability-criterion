## Results for Section 2.7.3: Cluster Analysis of Single-cell RNA Sequence Data
## Date: 4/20/2024

library(dplyr)
library(Seurat)

rm(list=ls())
source("code/PHM_algorithm.R")
source("code/visualizations.R")

pbmc_data <- Read10X(data.dir = "/mnt/c/Users/aturf/Downloads/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)

## Quality control steps
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## Normalize / scale all variables
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))


## Principal Component Analysis & Select number of PCs
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

screeplot <- ElbowPlot(pbmc, reduction="pca") +
  geom_point(aes(x=dims, y=stdev)) + geom_line(aes(x=dims, y=stdev)) +
  ylab("Standard Deviation") + xlab("") +
  scale_x_continuous(breaks=(0:10)*2) +
  theme_bw() + theme(panel.grid.minor=element_blank())

plot(screeplot)

ggsave("plots/pbmc_screeplot.png",
       plot=screeplot,
       height=3, width=6, units="in")


## Run PHM Procedure based on 10 PCs
num_pcs <- 10
embed_mat <- Embeddings(pbmc, reduction = "pca")
rownames(embed_mat) <- NULL
colnames(embed_mat) <- NULL
embed_mat <- embed_mat[, 1:num_pcs]

set.seed(20240308)
gmm_pbmc <- Mclust(embed_mat, G=1:15, warn=T)
phm_res <- PHM(gmm_pbmc, mc.samples=1e6, num.cores=7)

dendro <- plotPHMDendrogram(phm_res)
plot(dendro)

ggsave("plots/pbmc_phm_dendro.png",
       plot=dendro,
       height=4, width=6, units="in")


## Use the Louvain algorithm to identify the cell types present in the sample
## TODO: Add distruct plot

