# MASLD Subgroup Analysis

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

This repository contains the R code and documentation for the analysis of Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) using data from the Mayo Clinic Biobank (MCB) and Tapestry Study cohorts. The study is structured into four principal sections: 1. Data preparation, 2. MASLD subgroup identification, 3. Subgroup-specific downstream analyses, and 4. validation of subgroups in new cohorts. 

> **Note:** Data from the MCB and Tapestry cohorts are **not publicly available** due to privacy restrictions.

## 🔬 Overview

This repository contains the complete analytical framework for identifying clinically distinct MASLD subgroups through latent class analysis. The study integrates clinical, genomic, and longitudinal data to characterize disease heterogeneity and validate findings across independent cohorts.

### Key Features

- **Latent Class Analysis**: Advanced statistical methods for subgroup identification
- **Multi-cohort Validation**: Cross-validation using MCB and Tapestry Study data
- **Downstream Analyses**: Comprehensive evaluation of clinical outcomes, genomic risk, and medication responses
- **Reproducibility**: Multiple membership assignment methods for new patient classification

## 🧬 Study Design

The analysis is structured into four principal components:

1. **Data Preparation**: Preprocessing and transformation of clinical variables
2. **Subgroup Identification**: Latent class analysis to identify distinct MASLD phenotypes
3. **Downstream Analyses**: Subgroup-specific evaluation of:
   - Clinical outcomes and disease progression
   - Polygenic risk scores (PRS)
   - Medication utilization patterns
   - Genomic associations
4. **Validation**: Independent cohort validation and membership assignment methods

## 💻 Installation

### Prerequisites

- R version ≥ 4.0 (tested on R 4.2.2)
- RStudio (recommended)

### Required R Packages

Install all required packages by running:

```r
# Core packages for data manipulation and visualization
pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")

# Latent class analysis and clustering
pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")

# Statistical modeling and survival analysis
pkgs_additional <- c("olsrr", "dplyr", "tidyr", "purrr", "tidyfit", "table1", 
                     "RcppArmadillo", "BranchGLM", "MASS", "tidycmprsk", 
                     "survival", "survminer", "scatterplot3d", "mlr3misc")

# Combine and install all packages
all_pkgs <- unique(c(pkgs_overall, pkgs_lca, pkgs_additional))
install.packages(all_pkgs)

# Load packages
sapply(all_pkgs, require, character.only = TRUE, quietly = TRUE)
```

## 📁 Repository Structure

```
.
├── data/                           # Data directory (not included - see Data Availability)
│   ├── raw/                        # Raw datasets
│   └── processed/                  # Processed datasets
├── codes/                          # Analysis scripts
│   ├── preprocessing.R             # Data preparation and variable transformation
│   ├── lca_subgroup_identification.R  # Latent class analysis
│   ├── membership_methods/         # Patient assignment algorithms
│   │   ├── probability_based_assignment.R
│   │   ├── centroid_based_assignment.R
│   │   └── core_points_based_assignment.R
│   ├── downstream_analyses/        # Subgroup-specific analyses
│   │   ├── clinical_outcome.R      # Clinical outcomes analysis
│   │   ├── disease_progression.R   # Longitudinal disease tracking
│   │   ├── genotyping.R            # Genomic associations
│   │   ├── medication.R            # Medication usage patterns
│   │   └── mre_analysis.R          # MR elastography analysis
│   ├── prs_analysis.R              # Polygenic risk score analysis
│   └── visualize_results.R         # Figure generation
└── results/                        # Output directory for results and figures
```

## 🔄 Analysis Pipeline

### 1. Data Preprocessing (`preprocessing.R`)

Prepares the dataset for latent class analysis:
- Converts continuous variables to categorical using clinically relevant cutoffs
- Handles missing data
- Creates indicator variables for LCA

### 2. Subgroup Identification (`lca_subgroup_identification.R`)

Implements latent class analysis:
- Identifies distinct MASLD subgroups in the development cohort
- Applies model to full MCB and Tapestry datasets
- Examines subgroup-specific distributions of comorbidities and clinical variables

### 3. Membership Assignment Methods (`membership_methods/`)

Three approaches for assigning new patients to subgroups:

- **Probability-based**: Uses posterior probabilities from class-conditional response distributions
- **Centroid-based**: Assigns patients based on distance to subgroup centroids
- **Core points-based**: Utilizes representative core points from each subgroup

### 4. Polygenic Risk Score Analysis (`prs_analysis.R`)

Analyzes PRS distributions across subgroups using different SNP sets in both development and validation cohorts.

### 5. Downstream Analyses

Comprehensive subgroup characterization:

| Analysis | File | Description |
|----------|------|-------------|
| Clinical Outcomes | `clinical_outcome.R` | Evaluates disease severity and comorbidity patterns |
| Disease Progression | `disease_progression.R` | Longitudinal tracking of liver disease progression |
| Genomic Analysis | `genotyping.R` | Genetic variant associations with subgroups |
| Medication Patterns | `medication.R` | Analyzes medication utilization by subgroup |
| MRE Analysis | `mre_analysis.R` | Magnetic resonance elastography findings |

### 6. Visualization (`visualize_results.R`)

Generates publication-quality figures:
- Box plots comparing clinical variables across subgroups
- Bar plots showing disease prevalence
- Statistical comparisons for MCB development, validation, and Tapestry cohorts



