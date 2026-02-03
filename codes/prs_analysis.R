# =============================================
# PRS Analysis
# =============================================
# This script analyzes the distribution of polygenic risk scores (PRS) based on different SNP sets across subgroups 
# in both development and validation cohorts.
#
# Objectives:
# - Evaluate PRS proportions in the MCB development cohort and validation cohort
# - Evaluate PRS proportions in the Tapestry dataset with three assignment methods
#
# =============================================
# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "stringr",
  "data.table", "reshape2", "ggrepel", "scales",
  "paletteer", "poLCA", "networkD3", "scatterpie",
  "corrplot", "tidyLPA", "table1", "caret",
  "randomForest", "Boruta", "mlbench", "caTools"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------- Load Data -------------------
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
    prs_final <- read.delim(file, sep = "\t", header = TRUE)
    data_found <- TRUE
    break
  }
}

case_all_final <- merge(
  merged_data,
  prs_final,
  by = "PATIENT_ID",
  all.x = TRUE,
  all.y = FALSE
)

# ------------------- Categorize PRS Scores -------------------
prs_cut <- function(x) {
  case_when(
    x >= quantile(x, 0.9, na.rm = TRUE) ~ "High",
    x >= quantile(x, 0.4, na.rm = TRUE) ~ "Intermediate",
    TRUE ~ "Low"
  )
}

case_all_final <- case_all_final %>%
  mutate(
    PRS_Obs_C97   = prs_cut(PRS_97_SNP_SCORE1_SUM),
    PRS_t2D_C46   = prs_cut(PRS_46_SNP_SCORE1_SUM),
    PRS_cvd_C330  = prs_cut(PRS_330_SNP_SCORE1_SUM),
    PRS_depre_C8  = prs_cut(PRS_8_SNP_SCORE1_SUM),
    PRS_NAFLD_C15 = prs_cut(PRS_15),
    PRS_NAFLD_C5  = prs_cut(PNPLA3_TM6SF2_MBOAT7_GCKR_HSD17B13),
    PRS_NAFLD_DISCORD = prs_cut(DISCORD_3_SNP_SUM),
    PRS_NAFLD_CORD     = prs_cut(CORD_2_SNP_SUM)
  )

# ------------------- Prepare Data for Plotting -------------------
gene_columns <- c("PRS_NAFLD_CORD", "PRS_NAFLD_DISCORD")
long_data <- case_all_final %>%
  pivot_longer(
    cols = gene_columns,
    names_to = "PRS",
    values_to = "Value"
  )

# ------------------- Compute Proportions per Subgroup -------------------
percentage_data <- long_data %>%
  filter(!is.na(Subgroup)) %>%
  group_by(PRS, Subgroup, Value) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(PRS, Subgroup) %>%
  mutate(
    total_count = sum(Count),
    Percentage = (Count / sum(Count)) * 100
  )

# ------------------- Create Table 1 Summary -------------------
x1 <- table1(~.| Subgroup,
             data = case_all_final[, c(66,103:116)] %>% filter(!is.na(Subgroup)),
             overall = TRUE,
             topclass = "Rtable1-zebra")
tab1_df <- as.data.frame(x1)

write.table(tab1_df, "All_PRS_bio_new.csv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ------------------- Convert High Status to Binary -------------------
long_data <- long_data %>%
  mutate(Status = ifelse(Value == "High", 1, 0))

long_data$Subgroup <- factor(long_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

# ------------------- Summarize Proportions with SE -------------------
summary_data <- long_data %>%
  group_by(Subgroup, PRS) %>%
  summarise(
    Count = sum(Status),
    Total = n(),
    Proportion = Count / Total,
    SE = 1.96 * sqrt((Proportion * (1 - Proportion)) / Total),
    .groups = "drop"
  ) %>%
  mutate(Percentage = Proportion * 100)

# ------------------- Generate Bar Plot -------------------
disease_labels <- c(
  "PRS_NAFLD_DISCORD" = "MASLD PRS-Discordant",
  "PRS_NAFLD_CORD"   = "MASLD PRS-Concordant"
)

summary_data$PRS <- factor(summary_data$PRS, levels = names(disease_labels))
p_values <- long_data[, c(66,109:110)] %>%
  filter(Subgroup %in% all_subgroups) %>%
  group_by(PRS) %>%
  group_modify(~ {
    
    disease_data <- drop_na(.x)  # Remove NAs
    
    # Generate all pairwise comparisons of subgroups
    pairwise_results <- map_dfr(combn(all_subgroups, 2, simplify = FALSE), function(subgroup_pair) {
      
      # Subset data for this subgroup pair
      comparison_data <- disease_data %>%
        filter(Subgroup %in% subgroup_pair)
      print(comparison_data)
      # Ensure table has both groups (avoid empty levels)
      table_data <- table(comparison_data$Status, comparison_data$Subgroup)
      print(table_data)
      # Skip if table is incomplete
      if (ncol(table_data) < 2 || nrow(table_data) < 2) {
        return(tibble(Subgroup_1 = subgroup_pair[1], Subgroup_2 = subgroup_pair[2], p_value = NA))
      }
      
      # Apply Fisher or Chi-square test
      p_value <- tryCatch({
        if (min(table_data) < 5) {
          fisher.test(table_data, simulate.p.value = TRUE, B = 1e5)$p.value
        } else {
          chisq.test(table_data)$p.value
        }
      }, error = function(e) NA)  # Return NA if error occurs
      print(p_value)
      # Store results
      tibble(
        Subgroup_1 = subgroup_pair[1],
        Subgroup_2 = subgroup_pair[2],
        p_value = p_value
      )
    })
    
    return(pairwise_results)
  }) %>%
  ungroup() %>%
  group_by(PRS) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))  # Apply Bonferroni correction

disease_labels <- c(
  "PRS_Obs_C97" = "Obesity",
  "PRS_t2D_C46" = "T2D",
  "PRS_depre_C8" = "Depression",
  "PRS_cvd_C330" = "CVD",
  "PRS_NAFLD_DISCORD" = "MASLD PRS-Discordant",
  "PRS_NAFLD_CORD" = "MASLD PRS-Concordant",
  "PRS_NAFLD_C15" = "MASLD PRS-15",
  "PRS_NAFLD_C5" = "MASLD PRS-5"
)
disease_labels <- c(
  "PRS_NAFLD_DISCORD" = "MASLD PRS-Discordant",
  "PRS_NAFLD_CORD" = "MASLD PRS-Concordant"
)
#View(case_all_final)
summary_data$PRS <- factor(summary_data$PRS, levels = names(disease_labels))
summary_data
#View(summary_data)
#View(na.omit(summary_data))
p <- ggplot(na.omit(summary_data), aes(y = PRS, x = Percentage, fill = Subgroup)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, color = "gray40") +  # Gray border
  geom_errorbar(aes(xmin = Percentage - (SE * 100), xmax = Percentage + (SE * 100)),
                position = position_dodge(width = 0.9),  color = "black") +
  scale_fill_manual(values = c("#8B6914", "#B0C4DE", "#BBFFFF", "darkseagreen1", "lightsalmon", "plum1")) +
  scale_y_discrete(labels = disease_labels) +  # Relabel y-axis instead of x-axis
  labs(y= "", x = "Percentage (%)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 45, hjust = 1, color = "black", size = 16),  # Adjusted for vertical orientation
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    axis.text.x = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    legend.position = c(0.88,0.2),
    legend.text = element_text(color = "black", size = 15)
  )+
  guides(fill = guide_legend(reverse = TRUE))

# ------------------- Save Outputs -------------------
ggsave("bio_all_PRS_all_n_l_new_all.png", plot = p, width = 7, height = 7, dpi = 400)
ggsave("bio_all_PRS_all_n_l_new_all.pdf", plot = p, width = 7, height = 7, dpi = 400)
write.csv(summary_data, "PRS_bio_propotion_results_all.csv", row.names = FALSE)

write.csv(p_values, "p_value_PRS_bio_propotion_results_all.csv", row.names = FALSE)
