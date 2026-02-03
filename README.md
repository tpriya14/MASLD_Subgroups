# MASLD Subgroup Analysis

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

This repository contains the R code and documentation for the analysis of Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) using data from the Mayo Clinic Biobank (MCB) and Tapestry Study cohorts. The study is structured into four principal sections: 1. Data preparation, 2. MASLD subgroup identification, 3. Subgroup-specific downstream analyses, and 4. validation of subgroups in new cohorts. 

> **Note:** Data from the MCB and Tapestry cohorts are **not publicly available** due to privacy restrictions.

## 🔬 Overview

This repository contains the complete analytical framework for identifying clinically distinct MASLD subgroups through latent class analysis. The study integrates clinical, genomic, and longitudinal data to characterize disease heterogeneity and validate findings across independent cohorts.

## 🧬 Study Design

The analysis is structured into four principal components:

1. **Data Preparation**: Preprocessing and categorization of clinical variables according to the cutoff value for different biomarkers
2. **Subgroup Identification**: [Latent class analysis](https://pubmed.ncbi.nlm.nih.gov/37353681/) to identify distinct MASLD subgroups
3. **Downstream Analyses**: Subgroup-specific evaluation of:
   - Genetic variants analysis
   - Intrahepatic and extrahepatic complex disease risk analysis
   - Longitudinal biomarker analysis
   - Treatment outcome analysis
4. **Membership assignment methods**: To accurately map new patients to LCA-derived subgroups from the MCB development set, we benchmarked three membership assignment methods. 
   - Centroid-based: Computes Euclidean distances between patient feature vectors and subgroup centroids, assigning patients to the nearest subgroup.
   - Probability-based: Calculates [posterior membership probabilities](https://github.com/dlinzer/poLCA) using class-conditional response distributions from the development cohort, assigning patients to the subgroup with highest probability.
   - Core points-based: Identifies representative core members within each subgroup (inspired by [DBSCAN](https://github.com/mhahsler/dbscan?tab=readme-ov-file#ref-ester1996density)), then assigns new patients based on distance to the mean feature values of these core points.
5. **Validation**: Independent cohort validation with three membership assignment methods

## 💻 Installation

### Prerequisites

- R version ≥ 4.0 (tested on R 4.2.2)
- RStudio (recommended)

### Required R Packages

Install all required packages by running the below provided scripts.

## 📁 Repository Structure

```
.
├── data/                           # Data directory (not included)
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
│   │   └── mre_analysis.R          # Magnetic Resonance Elastography (MRE) analysis
│   ├── prs_analysis.R              # Polygenic risk score analysis
│   └── visualize_results.R         # Figure generation
└── results/                        # Output directory for results and figures
```

**References:**
Ester, M., Kriegel, H.-P., Sander, J. & Xu, X. A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise. Second International Conference on Knowledge Discovery and Data Mining (KDD'96), 226-231 (1996).

Zhou, W., et al. Latent class analysis-derived classification improves the cancer-specific death stratification of molecular subtyping in colorectal cancer. NPJ Precis Oncol 7, 60 (2023).
