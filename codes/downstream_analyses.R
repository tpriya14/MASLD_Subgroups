# Script: Downstream Analysis of Genotype Distribution, Disease Progression, and Medication Usage in MAYO Biobank Data
# Purpose: Perform downstream analysis of genotype distributions (PNPLA3, TM6SF2, HSD17B13, GCKR, MBOAT7),
#          disease progression (survival analysis for liver-related outcomes), and medication usage
#          across subgroups in the MAYO Biobank dataset, with statistical comparisons and visualizations.


# Load required libraries
library(ggsci)          # Color scales for scientific visualization
library(ggplot2)        # Data visualization
library(olsrr)          # Stepwise regression
library(scatterplot3d)  # 3D scatter plots
library(MASS)           # Statistical functions
library(dplyr)          # Data manipulation
library(tidyr)          # Data tidying
library(purrr)          # Functional programming
library(stringr)        # String manipulation
library(tidyfit)        # Model fitting
library(table1)         # Summary tables
library(caTools)        # Tools for data analysis
library(Boruta)         # Feature selection
library(mlbench)        # Machine learning benchmarks
library(caret)          # Machine learning utilities
library(randomForest)   # Random forest modeling
library(knitr)          # Dynamic report generation
library(survival)       # Survival analysis
library(tibble)         # Modern data frames
library(lubridate)      # Date handling
library(ggsurvfit)      # Survival plots
library(gtsummary)      # Summary tables for statistical models
library(tidycmprsk)     # Competing risks analysis
library(ggpattern)      # Patterned geoms for ggplot2
library(viridis)        # Color scales for heatmaps
library(cmprsk)         # Competing risks analysis
library(data.table)     # Efficient data handling
library(reshape2)       # Data reshaping
library(ggrepel)        # Avoid overlapping labels in plots
library(scales)         # Scale functions for visualization
library(paletteer)      # Comprehensive color palettes
library(poLCA)          # Latent Class Analysis
library(networkD3)      # Network visualization
library(scatterpie)     # Scatter pie plots
library(corrplot)       # Correlation plots
library(tidyLPA)        # Tidy Latent Profile Analysis
library(rstatix)        # Dunn's test and statistical annotations
library(ggpubr)         # Boxplots with statistical annotations

# Set working directory for output files
setwd("/MAYO_BIOBANK/Journal_Final_Pic/75")

# --- Part 1: Genotype Distribution Analysis ---

# Load dataset
merged_data <- read.delim("T100_5_Tap_After_KMeans_NAssign_matched_Bio_T2560_clinical_info_10_28_genetic.csv", 
                          sep = "\t", header = TRUE)

# Define gene columns for analysis
gene_columns <- c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7")

# Reshape data to long format for genotype analysis
long_data <- merged_data %>%
  pivot_longer(cols = all_of(gene_columns),
               names_to = "Gene",
               values_to = "Genotype")

# Calculate genotype percentages by gene and subgroup
percentage_data <- long_data %>%
  group_by(Gene, Subgroup, Genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Gene, Subgroup) %>%
  mutate(
    total_count = sum(Count),
    Percentage = (Count / total_count) * 100
  )

# Convert genotypes to binary status
long_data <- long_data %>%
  mutate(Status = ifelse(Genotype == "Yes", 1, 0))

# Ensure Subgroup is a factor with ordered levels
long_data$Subgroup <- factor(long_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

# Summarize genotype prevalence
summary_data <- long_data %>%
  group_by(Subgroup, Gene) %>%
  summarise(
    Count = sum(Status),
    Total = n(),
    Proportion = Count / Total,
    SE = 1.96 * (sqrt((Proportion * (1 - Proportion)) / Total)),
    .groups = "drop"
  ) %>%
  mutate(Percentage = Proportion * 100)

# Calculate p-values for genotype prevalence
all_subgroups <- c("Control", "C1", "C2", "C3", "C4", "C5")
p_values <- long_data %>%
  filter(Subgroup %in% all_subgroups) %>%
  group_by(Gene) %>%
  group_modify(~ {
    disease_data <- drop_na(.x)
    pairwise_results <- map_dfr(combn(all_subgroups, 2, simplify = FALSE), function(subgroup_pair) {
      comparison_data <- disease_data %>%
        filter(Subgroup %in% subgroup_pair)
      table_data <- table(comparison_data$Status, comparison_data$Subgroup)
      if (ncol(table_data) < 2 || nrow(table_data) < 2) {
        return(tibble(Subgroup_1 = subgroup_pair[1], Subgroup_2 = subgroup_pair[2], p_value = NA))
      }
      p_value <- tryCatch({
        if (min(table_data) < 5) {
          fisher.test(table_data, simulate.p.value = TRUE, B = 1e5)$p.value
        } else {
          chisq.test(table_data)$p.value
        }
      }, error = function(e) NA)
      tibble(
        Subgroup_1 = subgroup_pair[1],
        Subgroup_2 = subgroup_pair[2],
        p_value = p_value
      )
    })
    return(pairwise_results)
  }) %>%
  ungroup() %>%
  group_by(Gene) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))

# Filter p-values for Control comparisons
p_values_genotype <- p_values %>%
  filter(Subgroup_1 == "Control")

# Create interaction term for Gene and Genotype
percentage_data$Gene_Genotype <- interaction(percentage_data$Gene, percentage_data$Genotype)

# Define ordered levels for heatmap
percentage_data$Gene_Genotype <- factor(percentage_data$Gene_Genotype, 
                                        levels = c("MBOAT7.TT (1/1)", "MBOAT7.CT (0/1)", "MBOAT7.CC (0/0)",
                                                   "GCKR.TT (1/1)", "GCKR.CT (0/1)", "GCKR.CC (0/0)",
                                                   "HSD17B13.TT (1/1)", "HSD17B13.AT (0/1)", "HSD17B13.AA (0/0)",
                                                   "TM6SF2.TT (1/1)", "TM6SF2.CT (0/1)", "TM6SF2.CC (0/0)",
                                                   "PNPLA3.GG (1/1)", "PNPLA3.CG (0/1)", "PNPLA3.CC (0/0)"))

# Filter for risk alleles
percentage_data <- percentage_data %>% 
  filter(Gene_Genotype %in% c("MBOAT7.TT (1/1)", "GCKR.TT (1/1)", 
                              "HSD17B13.TT (1/1)", "TM6SF2.TT (1/1)", 
                              "PNPLA3.GG (1/1)"))

percentage_data$Subgroup <- factor(percentage_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

# Create heatmap for genotype distribution
heatmap_plot <- ggplot(percentage_data, aes(x = Subgroup, y = Gene_Genotype, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage), 
                color = ifelse(Percentage > 15, "white", "black")), 
            size = 4) +
  scale_fill_viridis(name = "Percentage %", option = "magma", direction = -1, limits = c(0, 30)) +
  scale_color_identity() +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save heatmap
ggsave("genotype_heatmap_bio.png", plot = heatmap_plot, width = 8, height = 4, dpi = 400)
ggsave("genotype_heatmap_bio.pdf", plot = heatmap_plot, width = 8, height = 4, dpi = 400)

# Export genotype results
write.table(percentage_data, file = "genotype_heatmap_bio.csv", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
write.table(p_values, file = "genotype_heatmap_bio_p_value.csv", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# --- Part 2: Disease Progression Analysis ---



# --- Part 3: Medication Analysis ---

