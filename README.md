# ONDSA: Omics Networks Differential and Similarity Analysis

ONDSA (Omics Networks Differential and Similarity Analysis) is a framework designed to analyze structural differences and similarities across multiple clinical conditions using continuous omics data. The goal of ONDSA is to provide insights into disease mechanisms by identifying differential and similar network structures across groups. It leverages the Gaussian Graphical Model (GGM), a statistical approach for representing conditional dependencies among components, enabling comprehensive exploration of disease mechanisms using multi-omics data. ONDSA conducts rigorous statistical tests to obtain differential and similar omics network structures   across groups with false discovery rate (FDR) control at both steps. Through simulations and application to multi-omics datasets, ONDSA has shown its utility and computational efficiency in biomarker discovery and pathway analysis in cognitive aging.

Here, we provide an overview of ONDSA:

<img src="ONDSAworkflow.png" alt="Overview of ONDSA" height="600">

## Features

- **Network Analysis**: Represents conditional dependencies using GGMs for disease mechanism exploration.
- **Multi-group Comparison**: Tests for differences and similarities across multiple biological conditions.
- **FDR Control**: Controls the FDR for robust statistical inference.
- **Designed for Multi-omics Data**: Suitable for continuous omics data, such as proteomics, metabolomics, etc.

## Installation

ONDSA and its experiments are implemented in R. Make sure the following dependencies are installed: "FastGGM", "Rcpp", and "RcppParallel". These packages are required to estimate precision matrices and run ONDSA effectively. You can install ONDSA using the following steps:

```r
library(devtools)
library(remotes)

# Install Rcpp and RcppParallel packages
install.packages("Rcpp")
install.packages("RcppParallel")

# Install FastGGM from GitHub
remotes::install_github("wt2015-github/FastGGM")

# Install ONDSA package from GitHub
remotes::install_github("jiachenchen322/ONDSA")
```
## Input Data Format
ONDSA requires estimated precision matrices from multiple groups as input. Each group's network should be represented by a precision matrix, which captures the conditional dependencies among the variables (e.g., genes, proteins). The precision matrices can be estimated from raw omics data using the FastGGM package, which provides efficient tools for precision matrix calculation. FastGGM can be used to estimate group-specific precision matrices. You could load your own data for multiple groups (rows representing individuals and columns representing variables) and calculate the precision matrices. GGMs are useful for modeling omics networks based on continuous variables, assuming that the data can be transformed to approximate normality. Appropriate transformations may be needed if the data deviate significantly from normality. ONDSA identifies both differential and similar network structures across groups, providing insights into the omics relationships among variables across different conditions.

## Usage

`ONDSA(Omega, p, K, N, alpha, varname)`

## Arguments

- **Omega**: A list of estimated precision matrices, each representing the omics network for a specific group.
- **p**: An integer, the number of variables (nodes) in each network.
- **K**: An integer, the number of groups (or networks).
- **N**: A numeric vector containing the sample sizes for each group.
- **alpha**: A numeric value specifying the FDR control level for differential and similar structures at both steps.
- **varname**: A character vector of names for the variables (nodes) in the networks.

## Details

The `ONDSA()` function performs a two-step differential and similarity analysis on multiple omics networks based on GGMs. It identifies both differential omics network structures (edges that differ across groups) and similar structures (edges that are similar across groups). The procedure conducts rigorous statistical tests with FDR control at both steps. ONDSA requires that each group is measured on the same set of variables, ensuring that the number of variables remains the same across all groups.

## Return

The `ONDSA()` function returns a list containing two data frames:

- **differential_structures**: A data frame of edges that differ significantly across the groups.
- **similar_structures**: A data frame of edges that are similar across the groups.

## Dependencies

- `FastGGM`: For efficient calculation of Gaussian Graphical Models.
- `Rcpp` and `RcppParallel`: Provide C++ interfaces to accelerate computations, used in conjunction with FastGGM for precision matrix estimation.

Make sure to install these dependencies before using ONDSA.

## Example Omics Data

The `ThreeGroupRawOmicsData` dataset contains raw omics data for three groups and can be used to estimate group-specific precision matrices. You could load your own data for multiple groups (rows representing individuals and columns representing variables) and calculate the precision matrices. GGMs are useful for modeling omics networks based on continuous variables, assuming that the data can be transformed to approximate normality. You can load the example dataset using the following command:

```r
data("ThreeGroupRawOmicsData", package = "ONDSA")
```

## Getting Started
Here is an example of using ONDSA to estimate precision matrices from raw omics data and perform differential and similarity analysis on multi-group omics networks.
```r
# Load the package
library(ONDSA)

# Load example data
# This dataset contains raw omics data for three groups
data("ThreeGroupRawOmicsData", package = "ONDSA")

# Standardize the raw data for each group
Data_1Std <- apply(group1_rawdata, 2, standardize)
Data_2Std <- apply(group2_rawdata, 2, standardize)
Data_3Std <- apply(group3_rawdata, 2, standardize)

# Set parameters
alpha <- 0.05
p <- 500
K <- 3
N <- c(nrow(group1_rawdata), nrow(group2_rawdata), nrow(group3_rawdata))
varname <- colnames(group1_rawdata)

# Estimate group-specific precision matrices using FastGGM
n_total <- sum(N)
lambda_value <- sqrt(2 * log(p / sqrt(n_total)) / n_total)
fastggm_1 <- FastGGM_Parallel(Data_1Std, lambda=lambda_value)
fastggm_2 <- FastGGM_Parallel(Data_2Std, lambda=lambda_value)
fastggm_3 <- FastGGM_Parallel(Data_3Std, lambda=lambda_value)

# Create list of precision matrices
Omega <- list(fastggm_1$precision, fastggm_2$precision, fastggm_3$precision)

# Run ONDSA with the estimated precision matrices
result <- ONDSA(Omega, p, K, N, alpha, varname)

# View differential and similar structures
print(result$differential_structures)
print(result$similar_structures)
```

## Citation

If you use ONDSA in your research, please cite our related publication/software.

## Contact

If you have any questions or suggestions, please create an issue on GitHub or contact the author/maintainer Jiachen Chen (chenjc@bu.edu).


