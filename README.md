# MASLD Study Analysis

This repository contains the R code and documentation for the analysis of Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) using data from the Mayo Clinic Biobank (MCB) and Tapestry Study cohorts. The study is structured into four principal sections: 1. data preparation, 2. MASLD subgroup identification, 3. subgroup-specific downstream analyses, and 4. validation of subgroups in new cohorts. Data from the MCB and Tapestry Study cohorts are not publicly available due to privacy restrictions.


## Required R Packages

The following R packages are required to run the analyses:

```R
# Install required packages
pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")
pkgs_additional <- c("olsrr", "dplyr", "tidyr", "purrr", "tidyfit", "table1", 
                     "RcppArmadillo", "BranchGLM", "MASS", "tidycmprsk", 
                     "survival", "survminer", "scatterplot3d", "mlr3misc")

# Combine all packages
all_pkgs <- unique(c(pkgs_overall, pkgs_lca, pkgs_additional))

# Install packages
install.packages(all_pkgs)

# Load packages
sapply(all_pkgs, require, character.only = TRUE, quietly = TRUE)
```

## Repository Structure

- `data/`: Placeholder for raw and processed datasets (not included due to privacy restrictions).
- `scripts/`: R scripts for data preprocessing, subgroup identification, statistical analyses, and visualization.
  - `preprocessing.R`: Data cleaning and harmonization using `dplyr`, `tidyr`, and `data.table`.
  - `subgroup_identification.R`: Latent class analysis and clustering using `poLCA` and `tidyLPA`.
  - `downstream_analyses.R`: Longitudinal risk, genomic, and medication analyses using `tidycmprsk`, `survival`, `survminer`, and `BranchGLM`.
  - `validation.R`: Subgroup assignment and reproducibility testing.
  - `visualizations.R`: Code for generating figures using `ggplot2`, `networkD3`, `scatterpie`, `corrplot`, `ggrepel`, and `scatterplot3d`.
- `results/`: Output directory for analysis results and figures.
- `README.md`: This file.

## Usage

1. **Setup**:
   - Ensure R version 4.2.2 is installed.
   - Install required packages as listed above.
   - Place raw data in the `data/` directory (ensure compliance with data access restrictions).

2. **Run Analyses**:
   - Execute scripts in the `scripts/` directory in the following order:
     1. `preprocessing.R`
     2. `subgroup_identification.R`
     3. `downstream_analyses.R`
     4. `validation.R`
     5. `visualizations.R`

3. **Output**:
   - Results (e.g., HRs, ORs, p-values) and figures will be saved in the `results/` directory.
