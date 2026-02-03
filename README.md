# MASLD Subgroup Analysis

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A comprehensive analysis pipeline for identifying and characterizing subgroups of Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) using latent class analysis and clinical data from the Mayo Clinic Biobank (MCB) and Tapestry Study cohorts.

## 📋 Table of Contents

- [Overview](#overview)
- [Study Design](#study-design)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Analysis Pipeline](#analysis-pipeline)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [Contact](#contact)

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

## 🚀 Usage

### Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/masld-analysis.git
   cd masld-analysis
   ```

2. **Install dependencies**:
   ```r
   source("install_packages.R")
   ```

3. **Prepare your data**:
   - Place raw data files in `data/raw/`
   - Ensure data format matches the expected structure

4. **Run the analysis pipeline**:
   Execute scripts in the following order:

   ```r
   # Step 1: Data preprocessing
   source("codes/preprocessing.R")
   
   # Step 2: Subgroup identification
   source("codes/lca_subgroup_identification.R")
   
   # Step 3: Downstream analyses
   source("codes/downstream_analyses/clinical_outcome.R")
   source("codes/downstream_analyses/disease_progression.R")
   source("codes/downstream_analyses/genotyping.R")
   source("codes/downstream_analyses/medication.R")
   
   # Step 4: Visualization
   source("codes/visualize_results.R")
   ```

5. **View results**:
   - Figures and tables will be saved in `results/`

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

## 📊 Data Availability

**Important**: The datasets from the Mayo Clinic Biobank (MCB) and Tapestry Study cohorts are **not publicly available** due to privacy restrictions and patient confidentiality requirements.

Researchers interested in accessing these data should contact:
- Mayo Clinic Biobank: [contact information]
- Tapestry Study: [contact information]

Data access is subject to institutional review board approval and data use agreements.

## 📝 Citation

If you use this code or methodology in your research, please cite:

```bibtex
@article{masld2024,
  title={Identification and Characterization of MASLD Subgroups through Latent Class Analysis},
  author={[Authors]},
  journal={[Journal]},
  year={2024},
  volume={[Volume]},
  pages={[Pages]},
  doi={[DOI]}
}
```

## 👥 Contact

For questions, issues, or collaboration inquiries:

- **Primary Investigator**: [Name] - [email]
- **Lead Analyst**: [Name] - [email]
- **Issues**: Please use the [GitHub Issues](https://github.com/yourusername/masld-analysis/issues) page

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Mayo Clinic Biobank participants and staff
- Tapestry Study participants and research team
- [Funding sources]

---

**Last Updated**: February 2026
