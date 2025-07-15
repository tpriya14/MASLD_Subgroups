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
- `codes/`: R scripts for data preprocessing, subgroup identification, statistical analyses, and visualization.
  - `preprocessing.R`: Data preparation before clustering.
  - `subgroup_identification.R`: Implementation of latent class analysis and subgroup identification.
  - `downstream_analyses.R`: Subgroup specific longitudinal risk, genomic, and medication analyses.
  - `validation.R`: Subgroup assignment and reproducibility testing in an independent dataset.
- `results/`: Output directory for analysis results and figures.

## Description of code files

| File name                                                     | Description                                               |
|---------------------------------------------------------------|-----------------------------------------------------------|
| [01_preprocessing.R](https://github.com/tpriya14/MASLD_Subgroups/blob/main/codes/preprocessing.R)                  | Mayo Biobank (MCB) and Tapestry Study data preparation                                      |
| [02_subgroup_identification.R](https://github.com/tpriya14/MASLD_Subgroups/blob/main/codes/lca_subgroup_identification.R)             | Latent class analysis was conducted on the development cohort to identify distinct subgroups. This analysis was also applied to the full MCB and Tapestry datasets to get reference subgroups. Additionally, the subgroup-specific distributions of comorbidities and clinical variables were examined.          
| [03_downstream_analyses.R](https://github.com/tpriya14/MASLD_Subgroups/blob/main/codes/downstream_analyses.R)                | Downstream analysis of genotype distribution, disease progression, and medication usage                    |
| [04_validation.R](https://github.com/tpriya14/MASLD_Subgroups/blob/main/codes/validation.R)        | Subgroup membership assignment methods                       |

## Usage

1. **Setup**:
   - Ensure R version 4.2.2 is installed.
   - Install required packages as listed above.
   - Place raw data in the `data/` directory.

2. **Run Analyses**:
   - Execute scripts in the `scripts/` directory in the following order:
     1. `preprocessing.R`
     2. `lca_subgroup_identification.R`
     3. `validation.R`
     4. `downstream_analyses.R`

3. **Output**:
   - Results and figures will be saved in the `results/` directory.
