# Implementation of the Distinguishability Criterion

This directory contains the main implementations of the methods and presented in the paper. The primary functions from each file are described below.

- [`Pmc_implementation.R`](Pmc_implementation.R)
    - Implementation of the Distinguishability criterion and the necessary helper functions
    - `computePmc()`: Implements calculation of the distinguishability criterion based on cubature integration
    - `computeMonteCarloPmc()`: Implements the calculation of the distinguishability criterion based on Monte Carlo integration
- [`PHM_algorithm.R`](PHM_algorithm.R)
    - Implement the functions necessary to run the Pmc Hierarchical Merging (PHM) algorithm
    - `PHM()`: Function to run the PHM algorithm
    - `computeDeltaPmcMatrix()`: Compute the matrix of Pmc reduction values for the merges using cubature integration to evaluate 
    - `computeMonteCarloDeltaPmcMatrix()`: Compute the matrix of Pmc reduction values for the merges using MC integration to evaluate 
- [`visualizations.R`](visualizations.R)
    - Implement the functions necessary for the main visualizations of the PHM algorithm. All functions take the output of the `PHM()` procedure as an argument
    - `plotPHMDendrogram`: Visualize the PHM procedure as a dendrogram
    - `plotPHMDistruct`: Generate the distruct plot (for the posterior assignment probabilities) for a given number of clusters at any point int he PHM procedure
    - `plotPmcMatrix`: Visualize the matrix of Pmc reduction values as a heatmap

