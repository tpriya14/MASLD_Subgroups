# =============================================================================
# MASLD Subgroup Comparison: Literature vs. LCA-Derived Clusters
# =============================================================================
# Description:
#   Assigns patients to literature-defined MASLD subgroups using robust
#   Euclidean distance from published medians/IQRs with appropriate unit conversion, 
#   then computes centroid distances between those subgroups and 
#   the LCA-derived clusters from this study. 
#   Produces heatmap visualizations for two external comparators:
#     1. Liu et al. (UK Biobank)  — 5 subgroups (PMID: 40640848)
#     2. Raverdy et al. (ABOS Cohort)               — 3 subgroups (PMID: 39653777)
# =============================================================================


# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "stringr", "tibble", "lubridate",
  "data.table", "reshape2", "scales", "paletteer", "ggrepel",
  "poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA",
  "table1", "caret", "randomForest", "Boruta", "mlbench", "caTools",
  "survival", "ggsurvfit", "gtsummary", "tidycmprsk", "survminer",
  "olsrr", "scatterplot3d", "MASS", "tidyfit", "feather",
  "RColorBrewer", "cluster"
)

new_pkgs <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_packages, library, character.only = TRUE))

dir.create("figures", showWarnings = FALSE)

# ------------------- Load Data -------------------
final_all_bio1 <- read.delim("datasets/MASLD_75_with_clusters.csv", sep = "\t", header = TRUE)
final_all_bio2 <- read.delim("datasets/MASLD_25_test_with_clusters.csv", sep = "\t", header = TRUE)

final_all_bio2 <- final_all_bio2[, -c(66)]
colnames(final_all_bio2)[colnames(final_all_bio2) == "LatentClassCluster1"] <- "LatentClassCluster"

merged_data <- rbind(final_all_bio1, final_all_bio2)

# ------------------- Define Variables -------------------
vari_lca  <- c("Age", "AVG_BMI", "ALT", "AST", "TRIGLYCERIDE", "LDL", "A1C")
cont_vars <- c("Age", "AVG_BMI", "ALT", "AST", "TRIGLYCERIDE", "LDL", "A1C")

# =============================================
# Part 1: Liu et al. (UK Biobank) Comparison
# =============================================

# ------------------- Define Liu Subgroup Parameters -------------------
groups <- c("Dyslipidemia", "Younger", "Obesity", "Inflammatory", "Hepatotoxic")
# medians_mat: Matrix of median values for each subgroup from the UK Biobank cohort
# iqrs_mat: Matrix of interquartile ranges (IQRs) for each subgroup, as reported in Hong et al. 

# ------------------- Define Assignment Function -------------------
assign_subgroup_relative <- function(patient_row, medians_mat, iqrs_mat, groups) {
  patient_vec <- c(
    as.numeric(patient_row$Age), as.numeric(patient_row$AVG_BMI),
    as.numeric(patient_row$ALT), as.numeric(patient_row$AST),
    as.numeric(patient_row$TRIGLYCERIDE),
    as.numeric(patient_row$LDL), as.numeric(patient_row$A1C)
  )
  group_distances <- sapply(1:5, function(i) {
    local_median <- medians_mat[i, ]
    local_iqr    <- iqrs_mat[i, ]
    z_local <- (patient_vec - local_median) / local_iqr
    sqrt(sum(z_local^2, na.rm = TRUE))
  })
  return(groups[which.min(group_distances)])
}

# ------------------- Assign Patients to Liu Subgroups -------------------
# medians_mat and iqrs_mat collected from the Hong et. al. paper to use here
merged_data$LSubgroup <- apply(merged_data, 1, function(row) {
  assign_subgroup_relative(as.list(row), medians_mat, iqrs_mat, groups)
})

table1(~ Age + AVG_BMI + ALT + AST + LDL + HDL + A1C + TRIGLYCERIDE | LSubgroup, merged_data)

# ------------------- Scale and Compute Centroids -------------------
selected_data <- merged_data
selected_data[cont_vars] <- scale(selected_data[cont_vars])

lca_centroids  <- aggregate(selected_data[, vari_lca], by = list(selected_data$Subgroup),  FUN = mean)
lca_centroids1 <- aggregate(selected_data[, vari_lca], by = list(selected_data$LSubgroup), FUN = mean)

lca_centroid_matrix     <- as.matrix(lca_centroids[, -1])
species_centroid_matrix <- as.matrix(lca_centroids1[, -1])

# ------------------- Compute Distance Matrix -------------------
distance_matrix <- matrix(NA, nrow = nrow(species_centroid_matrix), ncol = nrow(lca_centroid_matrix))

for (i in 1:nrow(species_centroid_matrix)) {
  for (j in 1:nrow(lca_centroid_matrix)) {
    distance_matrix[i, j] <- dist(rbind(lca_centroid_matrix[j, ], species_centroid_matrix[i, ]))
  }
}

# ------------------- Format Distance Matrix -------------------
correlation_matrix <- distance_matrix
colnames(correlation_matrix) <- rev(c(
  "C1: Non-obese cardiometabolic",
  "C2: Male-predominant cardiorenal",
  "C3: Female-predominant with \nobesity and mood disorders",
  "C4: Polygenic MASLD",
  "C5: Polygenic MASH"
))
rownames(correlation_matrix) <- c("Dyslipidemia", "Hepatotoxic", "Inflammatory", "Obesity", "Younger")

# ------------------- Plot Heatmap -------------------
correlation_long <- melt(correlation_matrix)
colnames(correlation_long) <- c("Cluster1", "Cluster2", "Distance")

min_val <- min(correlation_long$Distance, na.rm = TRUE)
max_val <- max(correlation_long$Distance, na.rm = TRUE)

correlation_long$text_color <- ifelse(
  abs(correlation_long$Distance) > 1.4 & abs(correlation_long$Distance) < 2.7,
  "black", "white"
)

brbg_colors <- brewer.pal(11, "BrBG")

p <- ggplot(correlation_long, aes(x = Cluster1, y = rev(Cluster2), fill = Distance)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = brbg_colors,
    name    = "Distance",
    breaks  = c(min_val, max_val),
    labels  = c("\u2191 High similarity", "\u2193 Low similarity"),
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      label.position = "bottom",
      barwidth       = 10,
      barheight      = 1
    )
  ) +
  geom_text(aes(label = round(Distance, 2), color = text_color), size = 6) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "UK Biobank cohort", y = "Our subgroups") +
  theme(
    axis.text          = element_text(size = 16, colour = "black"),
    axis.text.x        = element_text(size = 16, colour = "black", angle = 90, hjust = 1, vjust = 1),
    axis.text.y        = element_text(size = 16, colour = "black"),
    axis.title         = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.x       = element_blank(),
    legend.position    = "top",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 14),
    legend.title.align = 0.5,
    legend.text.align  = 0.5
  )

print(p)

ggsave("figures/Liu_vs_our_change_z_new.png", plot = p, width = 8, height = 5, dpi = 500)
ggsave("figures/Liu_vs_our_change_z_new.pdf", plot = p, width = 8, height = 5, dpi = 500)

# Perfomormed the same analysis with the medians_mat and iqrs_mat collected from the Raverdy et. al. paper 
