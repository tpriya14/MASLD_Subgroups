##############################################################
# PRS Analysis Pipeline
# This script analyzes the distribution of polygenic risk scores (PRS) based on different SNP sets across subgroups 
# in both development and validation cohorts.
#
# Objectives:
# - Evaluate PRS proportions in the MCB development cohort and validation cohort
# - Evaluate PRS proportions in the Tapestry dataset with three assignment methods

##############################################################

# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "stringr",
  "data.table", "reshape2", "ggrepel", "scales",
  "paletteer", "poLCA", "networkD3", "scatterpie",
  "corrplot", "tidyLPA", "table1", "caret",
  "randomForest", "Boruta", "mlbench", "caTools"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------- Load Data -------------------
cat("Step 1: Loading data with cluster assignments\n")
cat("----------------------------------------------\n")

# List of possible dataset files to load
data_files <- c(
  "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv",
  "datasets/prob_validation_tapestry_with_predictions.csv",
  "datasets/tapestry_dbscan_assignments.csv",
  "datasets/centroid_mcb_test_centroid_assignments.csv"
)

data_found <- FALSE

# Try to load the first available data file
for(file in data_files) {
  if(file.exists(file)) {
    merged_data <- read.delim(file, sep = "\t", header = TRUE)
    prs_final <- read.delim(file, sep = "\t", header = TRUE)
    cat("✓ Loaded data from:", file, "\n")
    data_found <- TRUE
    break
  }
}

# Merge the two datasets by PATIENT_ID
case_all_final = merge(
  merged_data,
  prs_final,
  by.x = "PATIENT_ID",
  by.y = "PATIENT_ID",
  all.x = TRUE,
  all.y = FALSE
)

# Display number of rows after merge
nrow(case_all_final)

# ------------------------------------------------------------
# Categorize PRS scores into High / Intermediate / Low
# ------------------------------------------------------------

# PRS 97 SNP obesity score categorization
case_all_final = case_all_final %>% mutate(
  PRS_Obs_C97 = case_when(
    PRS_97_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_97_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    PRS_97_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_97_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# PRS for Type 2 Diabetes
case_all_final = case_all_final %>% mutate(
  PRS_t2D_C46 = case_when(
    PRS_46_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_46_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    PRS_46_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_46_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# PRS for cardiovascular disease
case_all_final = case_all_final %>% mutate(
  PRS_cvd_C330 = case_when(
    PRS_330_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_330_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    PRS_330_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_330_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# PRS for depression
case_all_final = case_all_final %>% mutate(
  PRS_depre_C8 = case_when(
    PRS_8_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_8_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    PRS_8_SNP_SCORE1_SUM >= quantile(case_all_final$PRS_8_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# PRS for MASLD (15 SNP model)
case_all_final = case_all_final %>% mutate(
  PRS_NAFLD_C15 = case_when(
    PRS_15 >= quantile(case_all_final$PRS_15, probs = 0.9, na.rm = TRUE) ~ "High",
    PRS_15 >= quantile(case_all_final$PRS_15, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# PRS for MASLD (5 SNP model)
case_all_final = case_all_final %>% mutate(
  PRS_NAFLD_C5 = case_when(
    PNPLA3_TM6SF2_MBOAT7_GCKR_HSD17B13 >= quantile(case_all_final$PNPLA3_TM6SF2_MBOAT7_GCKR_HSD17B13, probs = 0.9, na.rm = TRUE) ~ "High",
    PNPLA3_TM6SF2_MBOAT7_GCKR_HSD17B13 >= quantile(case_all_final$PNPLA3_TM6SF2_MBOAT7_GCKR_HSD17B13, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# Discordant PRS score
case_all_final = case_all_final %>% mutate(
  PRS_NAFLD_DISCORD = case_when(
    DISCORD_3_SNP_SUM >= quantile(case_all_final$DISCORD_3_SNP_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    DISCORD_3_SNP_SUM >= quantile(case_all_final$DISCORD_3_SNP_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# Concordant PRS score
case_all_final = case_all_final %>% mutate(
  PRS_NAFLD_CORD = case_when(
    CORD_2_SNP_SUM >= quantile(case_all_final$CORD_2_SNP_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
    CORD_2_SNP_SUM >= quantile(case_all_final$CORD_2_SNP_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
    .default = "Low"
  )
)

# Display column names for verification
colnames(case_all_final)

# ------------------------------------------------------------
# Prepare data for plotting
# ------------------------------------------------------------

# Select PRS columns of interest
gene_columns <- c("PRS_NAFLD_CORD", "PRS_NAFLD_DISCORD")

# Convert data to long format for plotting
long_data <- case_all_final %>%
  pivot_longer(
    cols = gene_columns,
    names_to = "PRS",
    values_to = "Value"
  )

# ------------------------------------------------------------
# Compute proportions per subgroup
# ------------------------------------------------------------
percentage_data <- long_data %>%
  filter(!is.na(Subgroup)) %>%
  group_by(PRS, Subgroup, Value) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(PRS, Subgroup) %>%
  mutate(
    total_count = sum(Count),
    Percentage = (Count / sum(Count)) * 100
  )

# ------------------------------------------------------------
# Create Table 1 summary
# ------------------------------------------------------------
x1 <- table1(~.| Subgroup,
             data = case_all_final[, c(66,103:116)] %>% filter(!is.na(Subgroup)),
             overall = TRUE,
             topclass = "Rtable1-zebra")

tab1_df <- as.data.frame(x1)

# Export Table 1 results
write.table(tab1_df, "All_PRS_bio_new.csv",
            sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)

# ------------------------------------------------------------
# Convert High status to binary (1/0)
# ------------------------------------------------------------
long_data <- long_data %>%
  mutate(Status = ifelse(Value == "High", 1, 0))

long_data$Subgroup <- factor(long_data$Subgroup,
                             levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

# ------------------------------------------------------------
# Summarize proportions with standard errors
# ------------------------------------------------------------
summary_data <- long_data %>%
  group_by(Subgroup, PRS) %>%
  summarise(
    Count = sum(Status),
    Total = n(),
    Proportion = Count / Total,
    SE = 1.96 * (sqrt((Proportion * (1 - Proportion)) / Total)),
    .groups = "drop"
  ) %>%
  mutate(Percentage = Proportion * 100)

# ------------------------------------------------------------
# Generate bar plot
# ------------------------------------------------------------
disease_labels <- c(
  "PRS_NAFLD_DISCORD" = "MASLD PRS-Discordant",
  "PRS_NAFLD_CORD" = "MASLD PRS-Concordant"
)

summary_data$PRS <- factor(summary_data$PRS, levels = names(disease_labels))

p <- ggplot(na.omit(summary_data),
            aes(y = PRS, x = Percentage, fill = Subgroup)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9),
           width = 0.8, color = "gray40") +
  geom_errorbar(
    aes(xmin = Percentage - (SE * 100),
        xmax = Percentage + (SE * 100)),
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  scale_fill_manual(values = c("#8B6914", "#B0C4DE", "#BBFFFF",
                               "darkseagreen1", "lightsalmon", "plum1")) +
  scale_y_discrete(labels = disease_labels) +
  labs(y= "", x = "Percentage (%)") +
  theme_minimal()

# Display plot
p

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
ggsave("bio_all_PRS_all_n_l_new_all.png", plot = p, width = 7, height = 7, dpi = 400)
ggsave("bio_all_PRS_all_n_l_new_all.pdf", plot = p, width = 7, height = 7, dpi = 400)

write.csv(summary_data, "PRS_bio_propotion_results_all.csv", row.names = FALSE)
write.csv(p_values, "p_value_PRS_bio_propotion_results_all.csv", row.names = FALSE)
