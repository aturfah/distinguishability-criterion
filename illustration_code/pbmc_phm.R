## Results for Section 2.7.3: Cluster Analysis of Single-cell RNA Sequence Data
## Date: 4/20/2024

library(dplyr)
library(Seurat)
library(gridExtra)

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

## Use the Louvain algorithm with bandwidth=0.5 to identify the cell types present in the sample
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

mapCellType <- function(vec) {
  sapply(vec, function(idx) cluster_cell_types[idx + 1])
}
cluster_cell_types <- c("Naive CD4+ T", 
                        "CD14+ Mono",
                        "Memory CD4+ T",
                        "B",
                        "CD8+ T",
                        "FCGR3A+ Mono",
                        "NK", 
                        "DC",
                        "Platelet")

cluster_levels <- c("Naive CD4+ T", "DC", "NK", "Memory CD4+ T",
                   "Platelet", "CD8+ T",
                   "B", "FCGR3A+ Mono", 
                   "CD14+ Mono")

raw_labels <- pbmc$seurat_clusters
mapped_labels <- mapCellType(as.numeric(as.character(raw_labels)))
mapped_labels <- factor(mapped_labels, levels=cluster_levels)


## Run PHM Procedure based on 10 PCs
num_pcs <- 10
embed_mat <- Embeddings(pbmc, reduction = "pca")
rownames(embed_mat) <- NULL
colnames(embed_mat) <- NULL
embed_mat <- embed_mat[, 1:num_pcs]

set.seed(20240308)
gmm_pbmc <- Mclust(embed_mat, G=1:15, warn=T)
phm_res <- PHM(gmm_pbmc, mc.samples=1e6, data=embed_mat, num.cores=7)

colors <- c("posterior.1"="royalblue3",
            "posterior.2"="olivedrab",
            "posterior.3"="magenta",
            "posterior.4"="coral", 
            "posterior.5"="gold",
            "posterior.6"="coral4",
            "posterior.7"="lightgray",
            "posterior.8"="plum",
            "posterior.9"="skyblue"
)

plt_dendro <- plotPHMDendrogram(phm_res, colors=colors)
plt_distruct <- plotPHMDistruct(phm_res, 
                                labels=mapped_labels, 
                                colors=colors,
                                include_title=T)

plt_joined <- arrangeGrob(plt_dendro, 
                          plt_distruct, 
                          nrow=2, ncol=1,
                          heights=list( unit(3, "in"), unit(1, "in")  ))

plot(plt_joined)

ggsave("plots/pbmc_phm.png",
       plot=plt_joined,
       width=6, height=4,
       units="in")


