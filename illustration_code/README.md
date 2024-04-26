# Scripts to generate results

Scripts are self-contained and can be run as-is. They should be run with the "working directory" set to the repository directory (i.e. the `code` folder should be visible). Larger-scale simulation output will be saved as `.RData` files to the `output/` directory, and loaded if available. Generated figures will be saved to the `plots/` directory. Tables are output as text.

## Main Body Figures

- Figure 1 (Section 2.4): [`baudry_gmm_example.R`](baudry_gmm_example.R)
- Figure 2 (Section 2.5): [`kmeans_example.R`](kmeans_example.R)
- Figure 3 (Section 2.6): [`hclust_hyp_test.R`](hclust_hyp_test.R) and [`hclust_power.R`](hclust_power.R)
- Figure 4 (Section 2.7.1): [`penguins_kmeans_hclust.R`](penguins_kmeans_hclust.R)
- Figure 5 (Section 2.7.2): [`hgdp_phm.R`](hgdp_phm.R)
- Figure 6 (Section 2.7.3): [`pbmc_phm.R`](pbmc_phm.R)

## Supplementary Figures

- Supplementary Figure 1 (Section 2.1): [`decision_rules.R`](decision_rules.R)
- Supplementary Figure 2 (Section 2.1): [`pmc_2gaus_visualize.R`](pmc_2gaus_visualize.R)
- Supplementary Figure 3 (Section 2.7.1): [`penguins_kmeans_hclust.R`](penguins_kmeans_hclust.R)
- Supplementary Figure 4 (Section 2.7.2): [`hgdp_phm.R`](hgdp_phm.R)
- Supplementary Figure 5 (Section 2.7.3): [`pbmc_phm.R`](pbmc_phm.R)

## Supplementary Tables

- Supplementary Table 1 (Section 2.5): [`kmeans_example.R`](kmeans_example.R)
- Supplementary Table 2 (Section 2.7.1): [`penguins_kmeans_hclust.R`](penguins_kmeans_hclust.R)
- Supplementary Table 3 (Section 2.7.1): [`penguins_kmeans_hclust.R`](penguins_kmeans_hclust.R)
- Supplementary Table 4 (Section 4.5): [`montecarlo_comparison.R`](montecarlo_comparison.R)

## Additional Files

Files that contain functions used across multiple analyses.

- [`hyp_test_helpers.R`](hyp_test_helpers.R): Functions to perform the Gao et al. and Pmc hierarchical clustering procedures, as well as the function to perform hierarchical clustering.
- [`stability_helpers.R`](stability_helpers.R): Implements the functions to calculate Pmc for a range of *K* values for k-means and hierarchical clustering. Also compute the Silhouette score, Stability based on ARI, and Prediction strength.
- [`prediction_strength_function.R`](prediction_strength_function.R): Slightly modified version of `fpc::prediction.strength()` to work with the hierarchical clustering and k-means functions used elsewhere in this repo.
