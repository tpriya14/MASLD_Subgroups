# =============================================
# Genetic Variant Heatmap by LCA Cluster
# =============================================
# This script creates a heatmap showing the distribution of
# genetic risk variants across latent class clusters.
# =============================================
# Objective:
# - To evaluate how key MASLD-associated genetic variants (PNPLA3, TM6SF2, HSD17B13, GCKR, MBOAT7) are 
#   distributed across latent subgroups.
# - To provide a visual comparison of genotype frequencies 
#   between Control and latent class groups (C1–C5).
# =============================================

# ------------------- Load Required Packages -------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "viridis")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------- Load Data -------------------
# Step 1: Loading data with cluster assignments

data_files <- c(
  "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv",
  "datasets/prob_validation_tapestry_with_predictions.csv",
  "datasets/tapestry_dbscan_assignments.csv",
  "datasets/centroid_mcb_test_centroid_assignments.csv"
)

data_found <- FALSE
for(file in data_files) {
  if(file.exists(file)) {
    merged_data <- read.delim(file, sep = "\t", header = TRUE)
    data_found <- TRUE
    break
  }
}

# ------------------- Prepare Data -------------------
# Ensure genetic columns are factors
gene_columns <- c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7")

for(gene in gene_columns) {
  merged_data[[gene]] <- as.factor(merged_data[[gene]])
}

# Ensure Subgroup is factor with correct levels
merged_data$Subgroup <- factor(
  merged_data$Subgroup,
  levels = c("Control", "C1", "C2", "C3", "C4", "C5")
)

# ------------------- Create Heatmap -------------------
# Step 2: Creating heatmap visualization

merged_data <- as.data.frame(merged_data)
long_data <- merged_data %>%
  pivot_longer(
    cols = gene_columns,
    names_to = "Gene",
    values_to = "Genotype"
  )

# Calculate percentages
percentage_data <- long_data %>%
  group_by(Gene, Subgroup, Genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Gene, Subgroup) %>%
  mutate(
    total_count = sum(Count),
    Percentage = (Count / sum(Count)) * 100
  )

# Create interaction variable and filter top risk genotypes
percentage_data$Gene_Genotype <- interaction(percentage_data$Gene, percentage_data$Genotype)
percentage_data$Gene_Genotype <- factor(percentage_data$Gene_Genotype, levels = c(
  'MBOAT7.TT (1/1)', 'MBOAT7.CT (0/1)', 'MBOAT7.CC (0/0)',
  'GCKR.TT (1/1)', 'GCKR.CT (0/1)', 'GCKR.CC (0/0)',
  'HSD17B13.TT (1/1)', 'HSD17B13.AT (0/1)', 'HSD17B13.AA (0/0)',
  'TM6SF2.TT (1/1)', 'TM6SF2.CT (0/1)', 'TM6SF2.CC (0/0)',
  'PNPLA3.GG (1/1)', 'PNPLA3.CG (0/1)', 'PNPLA3.CC (0/0)'
))
percentage_data <- percentage_data %>% filter(
  Gene_Genotype %in% c('MBOAT7.TT (1/1)', 'GCKR.TT (1/1)', 'HSD17B13.TT (1/1)', 'TM6SF2.TT (1/1)', 'PNPLA3.GG (1/1)')
)

percentage_data$Subgroup <- factor(percentage_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

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
    panel.border = element_rect(color = "black", fill = NA)
  )

# ------------------- Save Results -------------------
dir.create("figures", showWarnings = FALSE)
dir.create("figures/genetics", showWarnings = FALSE)

ggsave("figures/genetics/genotype_heatmap_LCA.png", plot = heatmap_plot, width = 8, height = 5, dpi = 400)
ggsave("figures/genetics/genotype_heatmap_LCA.pdf", plot = heatmap_plot, width = 8, height = 5)

write.table(percentage_data, "figures/genetics/genotype_frequencies_by_cluster.csv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------- Statistical Analysis -------------------
# Chi-square tests for each gene across subgroups
stat_results <- data.frame()

for(gene in gene_columns) {
  gene_data <- merged_data[!is.na(merged_data$Subgroup) & !is.na(merged_data[[gene]]), ]
  cont_table <- table(gene_data$Subgroup, gene_data[[gene]])
  
  if(min(cont_table) >= 5) {
    test <- chisq.test(cont_table)
    test_name <- "Chi-square"
  } else {
    test <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    test_name <- "Fisher's exact"
  }
  
  stat_results <- rbind(stat_results, data.frame(
    Gene = gene,
    Test = test_name,
    p_value = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No")
  ))
}

write.csv(stat_results, "figures/genetics/genotype_statistical_tests.csv", row.names = FALSE)
print(heatmap_plot)

