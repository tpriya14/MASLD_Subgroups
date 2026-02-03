# =============================================
# Medication impact analysis
# =============================================
# This script analyzes medication usage patterns, their association 
# with disease outcome and impact on clinical biomarker across
# subgroups in the MCB an Tapestry dataset
#   dataset. 
#
# Objectives:
#   - Analyze medication usage patterns across patient subgroups
#   - Evaluate disease associations using logistic regression 
#     with adjusted and unadjusted models
#   - Generate forest plots for odds ratio visualization
#   - Perform longitudinal laboratory value analysis 
#     (HbA1C, ALT, BMI, HDL)
#   - Compare treatment responses across medication combinations
#   - Conduct statistical testing with Bonferroni correction
#
# =============================================
# Load required libraries
required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "stringr", "tibble", "lubridate",
  "ggplot2", "ggsci", "ggrepel", "ggpubr", "ggpattern", "viridis", "scales", "paletteer", "scatterpie", "networkD3", "corrplot", "scatterplot3d",
  "caret", "randomForest", "MASS", "olsrr", "Boruta", "mlbench", "tidyfit", "survival", "ggsurvfit", "cmprsk", "tidycmprsk", "poLCA", "tidyLPA", "rstatix", "table1"
)

# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "data.table","dplyr","tidyr","purrr","stringr","tibble","lubridate",
  "ggplot2","ggpubr","ggsci","ggrepel","ggpattern","viridis","scales","paletteer",
  "scatterpie","networkD3","corrplot","scatterplot3d","caret","randomForest","MASS",
  "olsrr","Boruta","mlbench","tidyfit","survival","ggsurvfit","cmprsk","tidycmprsk",
  "poLCA","tidyLPA","rstatix","table1"
)

new_pkgs <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_packages, library, character.only = TRUE))


##############################################################
# Part 1: Medication usage pattern
##############################################################

med_files <- list(
  list(file="Med_WEIGHT_LOSS_MED.csv", type="Weight Loss"),
  list(file="Med_STEROID.csv", type="Steroid"),
  list(file="Med_INSULIN.csv", type="Insulin"),
  list(file="Med_ANTIDIABETIC.csv", type="Antidiabetic"),
  list(file="Med_CHOLESTEROL_MED.csv", type="Cholesterol"),
  list(file="Med_ANTIDEPRESSANT.csv", type="Antidepressant"),
  list(file="Med_HYPERTENSION.csv", type="Hypertension")
)

all_med_data <- lapply(med_files, function(file) {
  read.delim(file, sep=",", header=TRUE)
}) %>% bind_rows()

# ------------------- Load MCB Patient Dataset -------------------
data_files <- c(
  "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv"
)

data_found <- FALSE
for(file in data_files) {
  if(file.exists(file)) {
    final_all_bio <- read.delim(file, sep = "\t", header = TRUE)
    cat("✓ Loaded patient data from:", file, "\n")
    data_found <- TRUE
    break
  }
}
if(!data_found) stop("No patient dataset found!")
merged_data <- merge(final_all_bio, all_med_data, by.x="PATIENT_ID", by.y="PATIENT_ID", all.x=TRUE)
total_patient_counts <- table(merged_data$Subgroup)

# Count percentages per Subgroup and MedName
med_summary <- merged_data %>%
  group_by(Subgroup, MedType, MedName) %>%
  summarise(Count = n_distinct(PATIENT_ID), .groups="drop") %>%
  mutate(Total = total_patient_counts[Subgroup],
         Percentage = (Count / Total) * 100)

# Create all combinations for heatmap
all_combinations <- expand.grid(
  Subgroup = unique(med_summary$Subgroup),
  Medication_Label = unique(med_summary$MedName),
  stringsAsFactors = FALSE
)

heatmap_data <- all_combinations %>%
  left_join(med_summary %>% select(Subgroup, Medication_Label = MedName, Percentage), 
            by = c("Subgroup", "Medication_Label")) %>%
  mutate(Percentage = ifelse(is.na(Percentage), 0, Percentage))

# Plot heatmap
heatmap <- ggplot(heatmap_data, aes(x=Subgroup, y=Medication_Label, fill=Percentage)) +
  geom_tile(color="white") +
  scale_fill_viridis(name="Usage %", option="magma", direction=-1, limits=c(0,80)) +
  geom_text(aes(label=round(Percentage,1), color=ifelse(Percentage>14,"white","black")), size=4) +
  scale_color_identity() +
  theme_minimal(base_size=14) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, size=12, face="bold"),
    axis.text.y = element_text(size=14),
    axis.title=element_blank(),
    legend.position="right",
    legend.title=element_text(size=14),
    legend.text=element_text(size=10)
  )

# Save plot
ggsave("medication_usage_heatmap.png", plot=heatmap, width=12, height=6, dpi=400)
ggsave("medication_usage_heatmap.pdf", plot=heatmap, width=12, height=6, dpi=400)

heatmap

# Save data
write.csv(all_med_data, "medication_usage_summary.csv", row.names = FALSE)


##############################################################
# Part 2: Disease associations with medication
##############################################################

factor_vars <- c("CVD1","CKD1","Fibrosis","Cirrhosis","Hepatocellular","Transplant",
                 "Depression1","Obesity1","Hyperlipidemia1","Diabetes1","Sleep1",
                 "Hypertension1","Joint1","Migraine1","Steatohepatitis")

merged_data[factor_vars] <- lapply(merged_data[factor_vars], as.factor)

# Composite liver disease
merged_data <- merged_data %>%
  mutate(ALD = as.factor(ifelse(Steatohepatitis=="1"|Hepatocellular=="1"|Fibrosis=="1"|Cirrhosis=="1", "1","0")))

merged_data$Subgroup <- factor(merged_data$Subgroup, levels=c('Control','C1','C2','C3','C4','C5'))
merged_data <- merged_data %>% filter(Subgroup != "Control")

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
write_csv(final_df, "results/all_forest_results.csv")

##############################################################
# Part 3: Medication impact on biomarker
##############################################################
run_lab_pipeline <- function(input_file, test_pattern, output_prefix) {
  
  followup_wt <- read.delim(input_file, sep=",", header=TRUE)
  
  alt_patients <- followup_wt %>%
    filter(str_detect(TEST_NAME, regex(test_pattern, ignore_case = TRUE)))
  
  alt_patients$Date_Time <- as.POSIXct(
    alt_patients$LAB_RESULT_DTM,
    format="%Y-%m-%d-%H.%M.%OS"
  )
  
  merged_data_wt_med <- merge(final_all_bio, alt_patients,
                              by.x = "PATIENT_ID", by.y = "PATIENT_ID", all.x = TRUE)
  
  merged_data_wt_med <- merged_data_wt_med %>%
    filter(year(LAB_RESULT_DTM) >= 2017 & year(LAB_RESULT_DTM) <= 2023)
  
  merged_data_wt_med <- merge(
    merged_data,
    merged_data_wt_med[, c(1,68:72)],
    by.x = "PATIENT_ID",
    by.y = "PATIENT_ID",
    all.x = TRUE
  )
  
  merged_data_wt_med$LAB_RESULT_DTM <- as.Date(merged_data_wt_med$LAB_RESULT_DTM)
  merged_data_wt_med$TEST_RESULT <- as.numeric(merged_data_wt_med$TEST_RESULT)
  
  data_last_10_years <- merged_data_wt_med
  
  med_combinations <- data_last_10_years %>%
    filter(Medication %in% c("Antidiabetic", "Insulin", "Cholesterol", "Weight Loss", "Corticosteroid")) %>%
    mutate(Medication = case_when(
      Medication %in% c("Antidiabetic", "Insulin") ~ "Antidiabetic/Insulin",
      TRUE ~ Medication
    )) %>%
    distinct(PATIENT_ID, Medication) %>%
    group_by(PATIENT_ID) %>%
    summarise(
      med_combo = paste(sort(unique(Medication)), collapse = " + ")
    ) %>%
    ungroup()

  
  data_with_med_combo <- data_last_10_years %>%
    left_join(med_combinations, by = "PATIENT_ID") %>%
    filter(med_combo %in% c(
      'Antidiabetic/Insulin + Cholesterol',
      'Antidiabetic/Insulin',
      'Cholesterol'
    ))
  
  date_diff_data <- data_with_med_combo %>%
    group_by(PATIENT_ID, med_combo, Subgroup) %>%
    summarise(
      first_value = TEST_RESULT[which.min(LAB_RESULT_DTM)],
      last_value  = TEST_RESULT[which.max(LAB_RESULT_DTM)],
      change = last_value - first_value,
      change_pct = (change / first_value) * 100
    ) %>%
    ungroup()

  
  date_diff_data2 <- date_diff_data %>%
    group_by(med_combo, Subgroup) %>%
    summarise(
      count = n(),
      first_mean = mean(first_value, na.rm = TRUE),
      last_mean  = mean(last_value, na.rm = TRUE),
      change_mean = mean(change, na.rm = TRUE),
      change_pct = mean(change_pct, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    group_by(Subgroup) %>%
    mutate(
      total_count = total_patient_counts[Subgroup],
      percentage = (count / total_count) * 100
    )

  
  stat.test <- date_diff_data %>%
    group_by(med_combo, Subgroup) %>%
    summarise(
      t_test_result = list(t.test(last_value, first_value, alternative = "two.sided", paired = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      p = sapply(t_test_result, function(x) x$p.value),
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

  
  date_diff_data2 <- date_diff_data2 %>%
    left_join(stat.test, by = c("Subgroup", "med_combo"))
  
  
  date_diff_data4 <- date_diff_data2 %>%
    select(-t_test_result)
  
  date_diff_data4$med_combo <- factor(
    date_diff_data4$med_combo,
    levels = c('Cholesterol', 'Antidiabetic/Insulin', 'Antidiabetic/Insulin + Cholesterol')
  )
  
  # FIXED: Removed duplicate left_join
  final_data <- date_diff_data4

  write.table(
    final_data,
    file = paste0(output_prefix, "_summary.csv"),
    sep = "\t",
    row.names = FALSE
  )
  

  p1 <- ggplot(date_diff_data4, aes(x = Subgroup, y = change_mean, fill = med_combo)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1), width = 0.8) +
    
    # Add change value labels
    geom_text(
      aes(label = ifelse(change_mean > 0, paste0('+', round(change_mean, 1)), round(change_mean, 1))), 
      position = position_dodge(width = 1),
      vjust = 0.5,
      hjust = ifelse(date_diff_data4$change_mean > 0, -0.1, 1),
      size = 4.5, 
      angle = 90
    ) +
    
    # Add baseline (first) value labels
    geom_text(
      aes(label = sprintf("%.1f", first_mean), y = 0), 
      position = position_dodge(width = 1),
      vjust = 0.5, 
      hjust = ifelse(date_diff_data4$change_mean > 0, 1.1, -0.1), 
      colour = "blue",
      size = 4.5,
      angle = 90
    ) +
    
    # Add count and percentage labels
    geom_text(
      aes(label = sprintf("%g (%.1f%%)", count, percentage), y = 0), 
      position = position_dodge(width = 1),
      vjust = 0.5, 
      hjust = ifelse(date_diff_data4$change_mean > 0, 1.7, -0.7), 
      colour = "red",
      size = 4,
      angle = 90
    ) +
    
    # Add significance labels inside bars
    geom_text(
      aes(label = ifelse(significance == "ns", '', significance), y = change_mean / 1.5),
      position = position_dodge(width = 1), 
      vjust = 0.5, 
      hjust = 0.5, 
      colour = "white", 
      size = 4.5
    ) +
    
    # Styling
    scale_shape_manual(values = c(21, 22, 23, 24)) + 
    labs(y = "Changes") +
    ylim(-4, 4) +
    scale_fill_jama() +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), color = "black", linetype = "dashed") +
    
    # Theme settings
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 20, color = "black"),
      legend.position = "top",
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black")
    )
  
  write.table(
    date_diff_data4, 
    file = paste0(output_prefix, "_plot_data.csv"), 
    sep = "\t", 
    col.names = TRUE, 
    quote = FALSE, 
    row.names = FALSE
  )
  
  ggsave(
    paste0(output_prefix, "_plot.png"), 
    plot = p1, 
    width = 8, 
    height = 4, 
    dpi = 500
  )
  
  return(final_data)
}
#########################################################
# RUN FOR ALL FOUR LABS
#########################################################

merged_data$ORDER_STOP_DTM  <- as.Date(merged_data$ORDER_STOP_DTM, format = "%Y-%m-%d")
merged_data$ORDER_START_DTM <- as.Date(merged_data$ORDER_START_DTM, format = "%Y-%m-%d")
merged_data$Date <- as.Date(merged_data$ORDER_START_DTM, format = "%Y-%m-%d")

setwd("follow_up/mcb/biobank_fold")
# A1C
result_A1C <- run_lab_pipeline(
  input_file = "A1C.csv",
  test_pattern = "HEMOGLOBIN|A1C|HGB",
  output_prefix = "A1C"
)

# ALT
result_ALT <- run_lab_pipeline(
  input_file = "ALT.csv",
  test_pattern = "ALT",
  output_prefix = "ALT"
)

# AST
result_AST <- run_lab_pipeline(
  input_file = "BMI.csv",
  test_pattern = "BMI",
  output_prefix = "BMI"
)

# HDL
result_HDL <- run_lab_pipeline(
  input_file = "HDL.csv",
  test_pattern = "HDL",
  output_prefix = "HDL"
)

#########################################################
# Repeat the same process with the Tapestry dataset
#########################################################


