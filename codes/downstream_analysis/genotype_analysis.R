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

# Load Required Packages
# ----------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "viridis")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

cat("========================================\n")
cat("Genetic Variant Heatmap Generation\n")
cat("========================================\n\n")

# ------------------- Load Data -------------------
cat("Step 1: Loading data with cluster assignments\n")
cat("----------------------------------------------\n")

# Try to load data with cluster assignments
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
    cat("✓ Loaded data from:", file, "\n")
    data_found <- TRUE
    break
  }
}


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

cat("✓ Data prepared\n")
cat("  Genes:", length(gene_columns), "\n")
cat("  Subgroups:", paste(levels(merged_data$Subgroup), collapse = ", "), "\n\n")



# ------------------- Create Heatmap -------------------
cat("Step 6: Creating heatmap visualization\n")
cat("---------------------------------------\n")

# Create the heatmap
merged_data$PNPLA3<- as.factor(merged_data$PNPLA3)
merged_data$GCKR <- as.factor(merged_data$GCKR)
merged_data$Subgroup <- as.factor(merged_data$Subgroup)
merged_data$TM6SF2 <- as.factor(merged_data$TM6SF2)
merged_data$HSD17B13 <- as.factor(merged_data$HSD17B13)
merged_data$MBOAT7 <- as.factor(merged_data$MBOAT7)
merged_data <- as.data.frame(merged_data)

gene_columns <- c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7")
long_data <- merged_data %>%
  pivot_longer(cols = c("PNPLA3",  "TM6SF2", "HSD17B13", "GCKR" , "MBOAT7"),
               names_to = "Gene", values_to = "Genotype")
# Calculate percentages
# Calculate percentages
percentage_data <- long_data %>%
  group_by(Gene, Subgroup, Genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Gene, Subgroup) %>%
  mutate(
    total_count = sum(Count),
    Percentage = (Count / sum(Count)) * 100) 


# Now create an interaction of Gene and Genotype to organize them accordingly
percentage_data$Gene_Genotype <- interaction(percentage_data$Gene, percentage_data$Genotype)


percentage_data$Gene_Genotype <- factor(percentage_data$Gene_Genotype, levels = c('MBOAT7.TT (1/1)', 'MBOAT7.CT (0/1)', 'MBOAT7.CC (0/0)','GCKR.TT (1/1)', 'GCKR.CT (0/1)', 'GCKR.CC (0/0)', 'HSD17B13.TT (1/1)', 'HSD17B13.AT (0/1)', 'HSD17B13.AA (0/0)',
                                                                                  'TM6SF2.TT (1/1)', 'TM6SF2.CT (0/1)', 'TM6SF2.CC (0/0)', 'PNPLA3.GG (1/1)', 'PNPLA3.CG (0/1)', 'PNPLA3.CC (0/0)'))
percentage_data <- percentage_data %>% filter(
  Gene_Genotype %in%  c('MBOAT7.TT (1/1)', 'GCKR.TT (1/1)', 'HSD17B13.TT (1/1)', 'TM6SF2.TT (1/1)', 'PNPLA3.GG (1/1)')
)
# Create the plot
percentage_data$Subgroup <- factor(percentage_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))
# Create the plot
heatmap_plot <- ggplot(percentage_data, aes(x = Subgroup, y = Gene_Genotype, fill = Percentage)) +
  geom_tile(color = "white") +  # Create the tiles with borders
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage), 
                color = ifelse(Percentage > 15, "white", "black")), 
            size = 4) +
  scale_fill_viridis(name = "Percentage %", option = "magma", direction = -1, limits = c(0, 30)) +
  scale_color_identity() + 
  #scale_fill_gradientn(colors = c("firebrick","darkgoldenrod1", "burlywood"), name = "Percentage") +  # Custom color gradient
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_rect(color = "black", fill = NA)  # Add border around the plot
  )

# Print the plot
print(heatmap_plot)

# Print the plot
print(heatmap_plot)
cat("✓ Heatmap created\n\n")

# ------------------- Save Results -------------------
cat("Step 7: Saving results\n")
cat("----------------------\n")

# Create output directory
dir.create("figures", showWarnings = FALSE)
dir.create("figures/genetics", showWarnings = FALSE)

# Save plots
ggsave(
  "figures/genetics/genotype_heatmap_LCA.png",
  plot = heatmap_plot,
  width = 8,
  height = 5,
  dpi = 400
)

ggsave(
  "figures/genetics/genotype_heatmap_LCA.pdf",
  plot = heatmap_plot,
  width = 8,
  height = 5
)

cat("✓ Heatmap saved:\n")
cat("  - figures/genetics/genotype_heatmap_LCA.png\n")
cat("  - figures/genetics/genotype_heatmap_LCA.pdf\n\n")

# Save data table
write.table(
  percentage_data,
  "figures/genetics/genotype_frequencies_by_cluster.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("✓ Data table saved:\n")
cat("  - figures/genetics/genotype_frequencies_by_cluster.csv\n\n")

# ------------------- Statistical Analysis -------------------
cat("Step 8: Performing statistical tests\n")
cat("-------------------------------------\n")

# Chi-square tests for each gene across subgroups
stat_results <- data.frame()

for(gene in gene_columns) {
  # Create contingency table
  gene_data <- merged_data[!is.na(merged_data$Subgroup) & 
                             !is.na(merged_data[[gene]]), ]
  
  cont_table <- table(gene_data$Subgroup, gene_data[[gene]])
  
  # Chi-square test
  if(min(cont_table) >= 5) {
    test <- chisq.test(cont_table)
    test_name <- "Chi-square"
  } else {
    test <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    test_name <- "Fisher's exact"
  }
  
  stat_results <- rbind(
    stat_results,
    data.frame(
      Gene = gene,
      Test = test_name,
      p_value = round(test$p.value, 4),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    )
  )
  
  cat(sprintf("  %s: p = %.4f (%s)\n", gene, test$p.value, test_name))
}

cat("\n✓ Statistical tests completed\n\n")

# Save statistical results
write.csv(
  stat_results,
  "figures/genetics/genotype_statistical_tests.csv",
  row.names = FALSE
)

cat("✓ Statistical results saved:\n")
cat("  - figures/genetics/genotype_statistical_tests.csv\n\n")

# ------------------- Summary Report -------------------
cat("========================================\n")
cat("Genetic Heatmap Complete!\n")
cat("========================================\n\n")

cat("Summary of genetic variant frequencies:\n")
cat("----------------------------------------\n\n")

for(gene in gene_columns) {
  cat(sprintf("%s:\n", gene))
  
  # Get risk allele frequency by subgroup
  gene_summary <- percentage_data %>%
    filter(grepl(gene, Gene_Genotype)) %>%
    select(Subgroup, Count, Percentage) %>%
    arrange(Subgroup)
  
  if(nrow(gene_summary) > 0) {
    for(i in 1:nrow(gene_summary)) {
      cat(sprintf("  %s: %d samples (%.1f%%)\n", 
                  gene_summary$Subgroup[i],
                  gene_summary$Count[i],
                  gene_summary$Percentage[i]))
    }
  }
  cat("\n")
}

cat("Key findings:\n")
cat("-------------\n")

# Identify highest risk genotype frequencies
high_risk <- percentage_data %>%
  group_by(Gene_Genotype) %>%
  slice_max(Percentage, n = 1) %>%
  ungroup()

cat("Highest risk allele frequencies:\n")
for(i in 1:nrow(high_risk)) {
  cat(sprintf("  %s: %s (%.1f%%)\n",
              high_risk$Gene_Genotype[i],
              high_risk$Subgroup[i],
              high_risk$Percentage[i]))
}

cat("\nSignificant associations (p < 0.05):\n")
sig_genes <- stat_results$Gene[stat_results$p_value < 0.05]
if(length(sig_genes) > 0) {
  for(gene in sig_genes) {
    p_val <- stat_results$p_value[stat_results$Gene == gene]
    cat(sprintf("  %s: p = %.4f\n", gene, p_val))
  }
} else {
  cat("  None\n")
}

cat("\n========================================\n")
cat("Files Generated:\n")
cat("========================================\n")


# Display the heatmap
print(heatmap_plot)
