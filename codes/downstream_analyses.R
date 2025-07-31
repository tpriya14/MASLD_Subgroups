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
merged_data <- read.delim("Tdataset.csv", 
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

# -------------------------------
# Define function (same as earlier)
# -------------------------------

run_survival_analysis <- function(
  merged_data,
  disease_file,
  disease_name,
  code_patterns,
  desc_keywords = NULL,
  output_prefix,
  date_column = "PROBLEM_NOTED_DTM.y",
  code_column = "dataset_CODE",
  desc_column = "DEPRESSION",
  subgroup_var = "Subgroup",
  sample_num_column = "sample_num",
  disease_sample_num_column = "sample_num"
) {
  disease_data <- read.delim(disease_file, sep = ",", header = TRUE)

  data <- merge(merged_data, disease_data, by.x = sample_num_column, by.y = disease_sample_num_column, all.x = TRUE)
  data$Date_d <- as.Date(data[[date_column]])

  # Dynamic disease status logic
  data <- data %>%
    mutate(
      DiseaseFlag = case_when(
        !is.null(desc_keywords) & grepl(paste(desc_keywords, collapse = "|"), .data[[desc_column]], ignore.case = TRUE) ~ "Yes",
        grepl(paste(code_patterns, collapse = "|"), .data[[code_column]]) ~ "Yes",
        .default = "No"
      ),
      Status = as.factor(ifelse(DiseaseFlag == "Yes", 1, 0))
    )

  followup_data <- data

  # Fit Cox model
  cox_formula <- as.formula(paste("Surv(Year1, Status) ~", subgroup_var))
  fit <- coxph(cox_formula, data = followup_data)
  cox_result <- tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = gsub(paste0(subgroup_var), "", term),
      HR = exp(estimate),
      lower_95 = exp(estimate - 1.96 * std.error),
      upper_95 = exp(estimate + 1.96 * std.error),
      p.value = sprintf("%.4f", p.value),
      HR_CI = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", lower_95), "-", sprintf("%.2f", upper_95), ")")
    ) %>%
    select(term, HR_CI, p.value)

  # Add baseline
  default_row <- data.frame(term = "Reference", HR_CI = "1.00", p.value = "-")
  summary_table <- rbind(default_row, cox_result)
  colnames(summary_table) <- c("Subgroup", "HR (95% CI)", "P Value")

  # Event summary
  event_summary <- followup_data %>%
    group_by(across(all_of(subgroup_var))) %>%
    summarise(
      Events = sum(Status == 1),
      Censored = sum(Status == 0),
      `% Events` = round(Events / (Events + Censored) * 100, 1),
      .groups = "drop"
    ) %>%
    rename(Subgroup = !!sym(subgroup_var))

  final_summary <- summary_table %>%
    left_join(event_summary, by = "Subgroup") %>%
    select(Subgroup, Events, Censored, `% Events`, `HR (95% CI)`, `P Value`)

  write.table(final_summary, paste0(output_prefix, "_summary.csv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Plot
  surv_plot <- cuminc(Surv(Year1, Status) ~ get(subgroup_var), data = followup_data) %>%
    ggcuminc(linewidth = 1) +
    labs(x = "Years") +
    ylim(0, 0.3) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    ) +
    add_risktable(
      risktable_stats = "n.risk",
      stats_label = list(n.risk = "Number at Risk"),
      size = 4.5,
      risktable_height = 0.25,
      theme = theme_risktable_default(axis.text.y.size = 12, plot.title.size = 14)
    )

  ggsave(paste0(output_prefix, "_survival_plot.png"), plot = surv_plot, width = 5, height = 6, dpi = 500)
  ggsave(paste0(output_prefix, "_survival_plot.pdf"), plot = surv_plot, width = 5, height = 6, dpi = 500)

  return(list(summary = final_summary, plot = surv_plot))
}

# -------------------------------
# Define ICD-10 Disease Code Map
# -------------------------------

disease_definitions <- list(
  "Sleep_Apnea" = list(codes = c("G47.3")),
  "Major_Depression" = list(codes = c("F32", "F33", "F34")),
  "Hypertension" = list(codes = c("I10")),
  "Fibrosis" = list(codes = c("K74.0")),
  "Cirrhosis" = list(codes = c("K74.6")),
  "HCC" = list(codes = c("C22")),
  "IHD" = list(codes = paste0("I", 20:25)),
  "ARF" = list(codes = c("N17")),
  "MASH" = list(codes = c("K75.8"))
)

# -------------------------------
# Run Analysis for All Diseases
# -------------------------------

for (disease in names(disease_definitions)) {
  cat("Running analysis for:", disease, "\n")
  
  run_survival_analysis(
    merged_data = merged_data,
    disease_file = "dataset_.csv",  # assumes all diseases are in this file
    disease_name = disease,
    code_patterns = disease_definitions[[disease]]$codes,
    desc_keywords = NULL,
    output_prefix = paste0("output/", disease)
  )
}

# --- Part 3: Longitudinal Biomarker Analysis: 2017 vs 2023 ---
# --------------------------------------------
# Biomarkers: BMI, ALT, AST, A1C, HDL
# --------------------------------------------

# Customize your file paths here
baseline_file <- "dataset.csv"
biomarker_file <- "-29-27.csv"  # Assumes all biomarkers are in the same long format
output_dir <- "output/"

# Load datasets
baseline <- read_csv(baseline_file)
biomarker_data <- read.delim(biomarker_file, sep = ",", header = TRUE)

# Convert date formats
biomarker_data$Date <- as.Date(biomarker_data$FLOWSHEET_ASSESSMENT_DTM)
biomarker_data$Year <- year(biomarker_data$Date)

# Merge with baseline demographics
merged <- merge(
  baseline,
  biomarker_data,
  by.x = "sample_num",
  by.y = "sample_num",
  all.x = FALSE
)

# Filter for relevant years and valid results
merged <- merged %>%
  filter(Year %in% c(2017, 2023)) %>%
  filter(FLOWSHEET_RESULT_TXT > 0 & FLOWSHEET_RESULT_TXT < 200)

# Biomarkers to analyze
biomarkers <- c("BMI", "ALT", "AST", "A1C", "HDL")

# Loop through each biomarker
for (bio in biomarkers) {
  message("Processing biomarker: ", bio)

  # Subset for current biomarker
  data_bio <- merged %>%
    filter(str_detect(FLOWSHEET_NAME, regex(bio, ignore_case = TRUE)))

  # Only include patients with data in both years
  common_ids <- data_bio %>%
    group_by(sample_num) %>%
    filter(n_distinct(Year) == 2) %>%
    pull(sample_num) %>%
    unique()

  data_bio <- data_bio %>% filter(sample_num %in% common_ids)

  # Calculate mean per subject per year
  summary_data <- data_bio %>%
    group_by(sample_num, Subgroup, Year) %>%
    summarise(mean_value = mean(FLOWSHEET_RESULT_TXT, na.rm = TRUE), .groups = "drop") %>%
    distinct()

  summary_data$Year <- as.factor(summary_data$Year)

  # Plot
  p <- ggerrorplot(
    summary_data, x = "Subgroup", y = "mean_value",
    desc_stat = "mean_sd", error.plot = "errorbar",
    color = "Year", merge = TRUE, add = "mean"
  ) +
    theme_minimal(base_size = 14) +
    labs(
      y = paste0("Mean ", bio),
      x = "Subgroup",
      title = paste0(bio, ": Change from 2017 to 2023")
    ) +
    theme(legend.position = "top")

  # Add lines for mean shifts
  mean_shift <- summary_data %>%
    group_by(Subgroup) %>%
    summarise(
      y2017 = mean(mean_value[Year == "2017"], na.rm = TRUE),
      y2023 = mean(mean_value[Year == "2023"], na.rm = TRUE),
      .groups = "drop"
    )

  p <- p + geom_segment(
    data = mean_shift,
    aes(
      x = as.numeric(Subgroup) - 0.2,
      xend = as.numeric(Subgroup) + 0.2,
      y = y2017,
      yend = y2023
    ),
    color = "blue"
  )

  # Paired t-test
  stat_test <- summary_data %>%
    group_by(Subgroup) %>%
    t_test(mean_value ~ Year, paired = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    mutate(
      p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
      p.adj = ifelse(p.adj < 0.001, "<0.001", sprintf("%.3f", p.adj)),
      Variable = bio
    )

  # Save results
  ggsave(paste0(output_dir, bio, ".png"), plot = p, width = 5, height = 4, dpi = 500)
  ggsave(paste0(output_dir, bio, ".pdf"), plot = p, width = 5, height = 4, dpi = 500)
  write.table(stat_test, paste0(output_dir, bio, "_ttest.tsv"), sep = "\t", row.names = FALSE)
}

# --- Part 3: Medication Analysis ---
# ---------------------------
# STEP 1: Prepare Data
# ---------------------------

# Convert medication columns to factor
med_columns <- c("Weight_Med", "Steroid_Med", "Insulin_Med", "Anti_Med", "Chol_Med", "AntiDep_Med", "Hyper_Med")
merged_data[med_columns] <- lapply(merged_data[med_columns], as.factor)

# Confirm structure
summary(merged_data)

# --------------------------------------------
# STEP 2: Define Helper Function to Calculate
#         Medication Use by Subgroup
# --------------------------------------------

calculate_med_usage <- function(df, med_var, med_type, total_counts) {
  df %>%
    group_by(Subgroup, !!sym(med_var)) %>%
    summarize(Unique_Patient_Count = n_distinct(sample_num), .groups = "drop") %>%
    mutate(
      Total_Count = total_counts[Subgroup],
      Percentage = (Unique_Patient_Count / Total_Count) * 100
    ) %>%
    ungroup() %>%
    mutate(MedType = med_type, MedName = !!sym(med_var)) %>%
    select(-!!sym(med_var))
}

# -------------------------------
# STEP 3: Total Patients per Subgroup
# (you must define this yourself)
# -------------------------------

# Example placeholder (replace with actual counts)
total_patient_counts <- merged_data %>%
  group_by(Subgroup) %>%
  summarize(count = n_distinct(sample_num)) %>%
  deframe()  # Converts to named vector

# -------------------------------
# STEP 4: Run Analysis per Med Type
# -------------------------------

result_list <- list(
  Weight = calculate_med_usage(merged_data, "Weight_Med", "Weight Loss", total_patient_counts),
  AntiDiabetic = calculate_med_usage(merged_data, "Anti_Med", "Antidiabetic", total_patient_counts),
  Insulin = calculate_med_usage(merged_data, "Insulin_Med", "Insulin", total_patient_counts),
  Steroid = calculate_med_usage(merged_data, "Steroid_Med", "Steroid", total_patient_counts),
  Cholesterol = calculate_med_usage(merged_data, "Chol_Med", "Cholesterol", total_patient_counts),
  Antidepressant = calculate_med_usage(merged_data, "AntiDep_Med", "Antidepressant", total_patient_counts),
  Hypertension = calculate_med_usage(merged_data, "Hyper_Med", "Hypertension", total_patient_counts)
)

# Combine all medication results
all_med_data <- bind_rows(result_list)

# -----------------------------------------
# STEP 5: Filter, Clean Labels, Prepare for Plotting
# -----------------------------------------

# Filter out missing or non-specific medication entries
all_med_data <- all_med_data %>%
  filter(!is.na(MedName) & MedName != "Others")

# Create label
all_med_data <- all_med_data %>%
  mutate(Medication_Label = as.character(MedName)) %>%
  mutate(Medication_Label = if_else(Medication_Label == "Insulin: Yes", "Insulin", Medication_Label))

# ----------------------------
# STEP 6: Prepare Data for Heatmap
# ----------------------------

# Ensure complete combinations of subgroup × medication
unique_subgroups <- unique(all_med_data$Subgroup)
unique_meds <- unique(all_med_data$Medication_Label)

all_combinations <- expand.grid(Subgroup = unique_subgroups, Medication_Label = unique_meds, stringsAsFactors = FALSE)

# Merge with actual results and fill missing with 0
heatmap_data <- all_combinations %>%
  left_join(all_med_data, by = c("Subgroup", "Medication_Label")) %>%
  mutate(Percentage = ifelse(is.na(Percentage), 0, Percentage))

# ----------------------------
# STEP 7: Create Heatmap
# ----------------------------

heatmap <- ggplot(heatmap_data, aes(x = Medication_Label, y = Subgroup, fill = Percentage)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Usage %", option = "magma", direction = -1, limits = c(0, 80)) +
  geom_text(
    aes(label = round(Percentage, 1),
        color = ifelse(Percentage > 14, "white", "black")),
    size = 4
  ) +
  scale_color_identity() +
  labs(x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )

# ----------------------------
# STEP 8: Save and View Output
# ----------------------------

ggsave("medication_usage_heatmap.png", heatmap, width = 10, height = 6, dpi = 400)
ggsave("medication_usage_heatmap.pdf", heatmap, width = 10, height = 6)

print(heatmap)

# Save data
write.csv(all_med_data, "medication_usage_summary.csv", row.names = FALSE)


#------------------------------------------
# 1. Preprocessing
#------------------------------------------

factor_vars <- c("CVD1", "CKD1", "Fibrosis", "Cirrhosis", "Hepatocellular", "Transplant",
                 "Depression1", "Obesity1", "Hyperlipidemia1", "Diabetes1", "Sleep1",
                 "Hypertension1", "Joint1", "Migraine1", "Steatohepatitis")

merged_data[factor_vars] <- lapply(merged_data[factor_vars], as.factor)

# Composite CLD variable
merged_data <- merged_data %>%
  mutate(CLD = case_when(
    Steatohepatitis == '1' ~ '1',
    Hepatocellular == '1' ~ '1',
    Fibrosis == '1' ~ '1',
    Cirrhosis == '1' ~ '1',
    TRUE ~ '0'
  )) %>%
  mutate(CLD = as.factor(CLD))


# Set subgroup as factor
merged_data$Subgroup <- factor(merged_data$Subgroup, levels = c('Control', 'C1', 'C2', 'C3', 'C4', 'C5'))
merged_data <- merged_data %>% filter(Subgroup != 'Control')

#------------------------------------------
# 2. Function to Create Forest Plot
#------------------------------------------

create_forest_plot <- function(data, disease, adjusted_model, unadjusted_model, medication) {
  # Model summaries
  adj_res <- tidy(adjusted_model, conf.int = TRUE) %>% filter(str_detect(term, "Subgroup"))
  unadj_res <- tidy(unadjusted_model, conf.int = TRUE) %>% filter(str_detect(term, "Subgroup"))

  # Clean term names
  adj_res$term <- gsub("Subgroup", "", adj_res$term)
  unadj_res$term <- gsub("Subgroup", "", unadj_res$term)

  # Bonferroni correction
  adj_res$p.adj <- p.adjust(adj_res$p.value, method = "bonferroni")
  unadj_res$p.adj <- p.adjust(unadj_res$p.value, method = "bonferroni")

  # Subgroup sizes
  counts <- data %>%
    group_by(Subgroup) %>%
    summarise(Total = n(), .groups = "drop")

  # Combine all results
  results <- tibble(
    Disease = disease,
    Medication = medication,
    Subgroup = c("C1", adj_res$term, adj_res$term),
    Model = c("Reference", rep("Unadjusted", nrow(unadj_res)), rep("Adjusted", nrow(adj_res))),
    OR = c(1, exp(unadj_res$estimate), exp(adj_res$estimate)),
    LowerCI = c(1, exp(unadj_res$conf.low), exp(adj_res$conf.low)),
    UpperCI = c(1, exp(unadj_res$conf.high), exp(adj_res$conf.high)),
    P_Value = c(NA, unadj_res$p.value, adj_res$p.value),
    P_Adj = c(NA, unadj_res$p.adj, adj_res$p.adj)
  ) %>%
    left_join(counts, by = "Subgroup") %>%
    mutate(Significance = ifelse(P_Adj < 0.05, "P < 0.05", "NS")) %>%
    drop_na(OR, LowerCI, UpperCI)

  # Plot
  forest_plot <- ggplot(results, aes(y = Subgroup, x = OR, xmin = LowerCI, xmax = UpperCI, color = Model)) +
    geom_pointrange(aes(size = Significance), position = position_dodge(width = 0.5)) +
    geom_errorbarh(position = position_dodge(width = 0.5), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_manual(values = c("Adjusted" = "red", "Unadjusted" = "blue")) +
    scale_size_manual(values = c("NS" = 0.5, "P < 0.05" = 1.5)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12, color = "black"),
      panel.grid = element_blank()
    ) +
    labs(x = "", y = "", color = "Model", size = "Significance")

  # Output paths
  out_dir <- file.path("results", "plots", medication)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  ggsave(filename = file.path(out_dir, paste0(medication, "_", disease, "_forest.png")),
         plot = forest_plot, width = 7, height = 6)
  ggsave(filename = file.path(out_dir, paste0(medication, "_", disease, "_forest.pdf")),
         plot = forest_plot, width = 7, height = 6)

  # Optional CSV
  # write_csv(results, file.path(out_dir, paste0(medication, "_", disease, "_results.csv")))

  return(results)
}

#------------------------------------------
# 3. Run Models for Multiple Diseases & Medications
#------------------------------------------

# Define diseases and medication predictors
diseases <- c("CLD", "CVD1", "CKD1", "Diabetes1", "Sleep1", "Depression1", "Transplant")
medications <- c("Weight_Loss", "Antidiabetic+Insulin+Cholestorol")

all_results <- list()

for (disease in diseases) {
  for (med in medications) {
    # Build models
    adj_model <- glm(as.formula(paste(disease, "~ Subgroup +", med)), data = merged_data, family = binomial)
    unadj_model <- glm(as.formula(paste(disease, "~ Subgroup")), data = merged_data, family = binomial)

    # Generate plot and results
    res <- create_forest_plot(merged_data, disease, adj_model, unadj_model, med)
    all_results[[paste(disease, med, sep = "_")]] <- res
  }
}

# Combine all results if needed
final_df <- bind_rows(all_results)
# write_csv(final_df, "results/all_forest_results.csv")


# --- Configuration ---
lab_keywords <- list(
  A1C = c("A1C"),
  ALT = c("ALT"),
  AST = c("AST"),
  HDL = c("HDL")
)

# --- Function to Process Lab Type ---
process_lab_measure <- function(lab_name, keywords) {

  # --- Load Lab Data ---
  lab_data <- read.delim(paste0("_", lab_name, "2024"), sep = ",", header = TRUE)
  filtered_lab <- lab_data %>%
    filter(str_detect(TEST_NAME, regex(paste(keywords, collapse = "|"), ignore_case = TRUE)))

  filtered_lab$Date_Time <- as.POSIXct(filtered_lab$LAB_RESULT_DTM, format = "%Y-%m-%d-%H.%M.%OS")
  filtered_lab$Month <- format(filtered_lab$Date_Time, "%Y-%m")

  merged <- merge(final_all_bio, filtered_lab, by.x = "sample_num", by.y = "sample_num", all.x = TRUE) %>%
    filter(year(LAB_RESULT_DTM) >= 2017 & year(LAB_RESULT_DTM) <= 2023)

  merged <- merge(merged_data, merged[, c(1, 68:72)], by = "sample_num", all.x = TRUE)
  merged$LAB_RESULT_DTM <- as.Date(merged$LAB_RESULT_DTM)
  merged$TEST_RESULT <- as.numeric(merged$TEST_RESULT)

  data_last_10_years <- merged

  # --- Medication Combinations ---
  med_combinations <- data_last_10_years %>%
    filter(Medication %in% c("Antidiabetic", "Insulin", "Cholesterol", "Weight Loss", "Corticosteroid")) %>%
    mutate(Medication = case_when(
      Medication %in% c("Antidiabetic", "Insulin") ~ "Antidiabetic/Insulin",
      TRUE ~ Medication
    )) %>%
    distinct(sample_num, Medication) %>%
    group_by(sample_num) %>%
    summarise(med_combo = paste(sort(unique(Medication)), collapse = " + ")) %>%
    ungroup() %>%
    add_count(med_combo, name = "num_patients")

  # --- Merge Med Combos and Filter ---
  data_with_med_combo <- data_last_10_years %>%
    left_join(med_combinations, by = "sample_num") %>%
    filter(med_combo %in% c('Antidiabetic/Insulin + Cholesterol', 'Antidiabetic/Insulin', 'Cholesterol'))

  # --- Calculate Changes ---
  date_diff_data <- data_with_med_combo %>%
    group_by(sample_num, med_combo, Subgroup) %>%
    summarise(
      min_date = min(LAB_RESULT_DTM, na.rm = TRUE),
      max_date = max(LAB_RESULT_DTM, na.rm = TRUE),
      med_start_date = as.Date(max(ifelse((min_date < ORDER_START_DTM & max_date > ORDER_START_DTM), ORDER_START_DTM, NA))),
      med_stop_date = as.Date(min(ifelse((year(ORDER_STOP_DTM) <= 2023 & ORDER_STOP_DTM > ORDER_START_DTM), ORDER_STOP_DTM, max_date))),
      Year1 = as.numeric(difftime(med_stop_date, med_start_date, units = "days")) / 365.25,
      first_val = TEST_RESULT[which.min(LAB_RESULT_DTM)],
      last_val = TEST_RESULT[which.max(LAB_RESULT_DTM)],
      change_pct = ((last_val - first_val) / first_val) * 100
    ) %>%
    ungroup() %>%
    group_by(sample_num) %>%
    filter(n_distinct(med_combo) == 1)

  # --- Summarize ---
  summary_data <- date_diff_data %>%
    group_by(med_combo, Subgroup) %>%
    summarise(
      count = n_distinct(sample_num),
      first_val_mean = mean(first_val),
      last_val_mean = mean(last_val),
      change = last_val_mean - first_val_mean,
      change_pct = ((last_val_mean - first_val_mean) / first_val_mean) * 100
    ) %>%
    ungroup() %>%
    group_by(Subgroup) %>%
    mutate(
      total_count = total_patient_counts[Subgroup],
      percentage = (count / total_count) * 100
    )

  # --- T-Test ---
  stat.test <- date_diff_data %>%
    group_by(med_combo, Subgroup) %>%
    summarise(t_test = list(t.test(last_val, first_val, paired = TRUE)), .groups = 'drop') %>%
    mutate(
      p = sapply(t_test, function(x) x$p.value),
      p.adj = p.adjust(p, method = "bonferroni"),
      p = case_when(
        p < 0.0001 ~ "<0.0001",
        p < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p)
      ),
      p.adj = case_when(
        p.adj < 0.0001 ~ "<0.0001",
        p.adj < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p.adj)
      ),
      significance = case_when(
        p == "<0.0001" ~ "****",
        p == "<0.001" ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  summary_data <- summary_data %>%
    left_join(stat.test, by = c("Subgroup", "med_combo")) %>%
    select(-t_test)

  summary_data$med_combo <- factor(summary_data$med_combo, levels = c('Cholesterol', 'Antidiabetic/Insulin', 'Antidiabetic/Insulin + Cholesterol'))

  # --- Plot ---
  plot <- ggplot(summary_data, aes(x = Subgroup, y = change, fill = med_combo)) +
    geom_bar(stat = "identity", position = position_dodge(1), width = 0.8) +
    geom_text(aes(label = ifelse(change > 0, paste0('+', round(change, 1)), round(change, 1))),
              position = position_dodge(1), vjust = 0.5, hjust = ifelse(summary_data$change > 0, -0.1, 1), size = 4.5, angle = 90) +
    geom_text(aes(label = sprintf("%.1f", first_val_mean), y = 0),
              position = position_dodge(1), hjust = 1.1, colour = "blue", size = 4.5, angle = 90) +
    geom_text(aes(label = sprintf("%g (%.1f%%)", count, percentage), y = 0),
              position = position_dodge(1), hjust = 1.7, colour = "red", size = 4, angle = 90) +
    geom_text(aes(label = ifelse(significance == "ns", '', significance), y = change / 1.5),
              position = position_dodge(1), vjust = 0.5, colour = "white", size = 4.5) +
    scale_fill_jama() +
    labs(y = paste0("Changes in ", lab_name)) +
    ylim(-4, 4) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dashed") +
    theme_minimal(base_size = 16) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    )

  ggsave(paste0("med_fp_cl/", lab_name, "_forest_plot.png"), plot = plot, width = 14, height = 5, dpi = 500)
  write.table(summary_data, file = paste0("med_fp_cl/", lab_name, "_summary.csv"), sep = ",", row.names = FALSE)
}

# --- Run for All Labs ---
setwd("/Biobank_fold")
for (lab in names(lab_keywords)) {
  process_lab_measure(lab, lab_keywords[[lab]])
}
