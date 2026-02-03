# =============================================
# Follow-Up Clinical Outcome Analysis
# =============================================
# This script processes longitudinal clinical measurements (BMI, ALT, HDL, A1C, AST, ALP)
# for patients in the MCB/Tapestry datasets.
#
# Objective:
# - Reads patient datasets and follow-up lab measurements
# - Calculates mean per patient per year
# - Creates paired plots showing the mean values with standard deviation (SD) 
#   and lines connecting paired observations.
# - Performs two sided paired t-tests with Bonferroni correction
# - Saves plots and statistics for all measures
=============================================

# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "ggplot2","ggpubr","dplyr","tidyr","purrr","stringr",
  "lubridate","rstatix","data.table","reshape2","ggrepel",
  "scales","paletteer","poLCA","networkD3","scatterpie",
  "corrplot","tidyLPA","table1","caret","randomForest",
  "Boruta","mlbench","caTools"
)

# Install missing packages and load all
new_pkgs <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_packages, library, character.only = TRUE))s

# ------------------- Analysis Function -------------------
# Function to process and plot paired measurements over years
analyze_followup <- function(data, measure_col, years = c(2017, 2023), 
                             subgroup_col = "Subgroup", file_prefix = "BMI",
                             ylabel = NULL) {
  
  # Convert date
  data$FLOWSHEET_ASSESSMENT_DTM <- as.Date(data$FLOWSHEET_ASSESSMENT_DTM)
  
  # Ensure measure is numeric
  data[[measure_col]] <- as.numeric(data[[measure_col]])
  
  # Filter for the selected years and valid values
  data_filtered <- data %>%
    filter(
      year(FLOWSHEET_ASSESSMENT_DTM) %in% years
    ) %>%
    mutate(
      year = year(FLOWSHEET_ASSESSMENT_DTM),
      !!sym(subgroup_col) := as.factor(!!sym(subgroup_col))
    )
  
  # Identify common IDs present in both years
  common_ids <- data_filtered %>%
    group_by(PATIENT_ID) %>%
    filter(n_distinct(year) == length(years)) %>%
    pull(PATIENT_ID) %>%
    unique()
  
  # Filter to only common IDs
  data_filtered <- data_filtered %>%
    filter(PATIENT_ID %in% common_ids)
  
  # Calculate mean per patient per year
  result_data <- data_filtered %>%
    group_by(PATIENT_ID, !!sym(subgroup_col), year) %>%
    summarise(mean_value = mean(!!sym(measure_col), na.rm = TRUE), .groups = "drop") %>%
    distinct() %>%
    mutate(year = as.factor(year))
  
  # Calculate mean values by subgroup for connecting lines
  mean_values <- result_data %>%
    group_by(!!sym(subgroup_col)) %>%
    summarise(
      mean_year1 = mean(mean_value[year == years[1]], na.rm = TRUE),
      mean_year2 = mean(mean_value[year == years[2]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create plot
  p1 <- ggerrorplot(
    result_data, 
    x = subgroup_col, 
    y = "mean_value", 
    desc_stat = "mean_sd",
    error.plot = "errorbar", 
    color = "year",
    merge = TRUE,
    add = "mean"
  ) +
    geom_segment(
      data = mean_values, 
      aes(
        x = as.numeric(!!sym(subgroup_col)) - 0.2, 
        xend = as.numeric(!!sym(subgroup_col)) + 0.2, 
        y = mean_year1, 
        yend = mean_year2
      ), 
      color = "blue"
    ) +
    labs(
      x = "Subgroups",
      y = if(is.null(ylabel)) paste0("Mean ", file_prefix) else ylabel,
      fill = "Year"
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "none"
    )
  
  # Statistical test
  stat_test <- result_data %>%
    group_by(!!sym(subgroup_col)) %>%
    t_test(mean_value ~ year, paired = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    left_join(mean_values, by = subgroup_col) %>%
    mutate(
      Variable = file_prefix,
      p = case_when(
        p < 0.0001 ~ "<0.0001",
        p < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p)
      ),
      p.adj = case_when(
        p.adj < 0.0001 ~ "<0.0001",
        p.adj < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p.adj)
      )
    )
  
  # Save statistical results
  write.table(
    stat_test, 
    paste0(file_prefix, "_stats.csv"), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
  )
  
  # Save plots
  ggsave(paste0(file_prefix, ".png"), plot = p1, width = 4, height = 4, dpi = 500)
  ggsave(paste0(file_prefix, ".pdf"), plot = p1, width = 4, height = 4, dpi = 500)
  
  return(list(plot = p1, stat_test = stat_test))
}


setwd("follow_up/mcb/Biobank_fold")

measure_config <- list(
  BMI = list(file = "BMI.csv", ylabel = expression(paste("Mean BMI (kg/m"^{2},")"))),
  ALT = list(file = "ALT.csv", ylabel = "Mean ALT (U/L)"),
  HDL = list(file = "HDL.csv", ylabel = "Mean HDL (mg/dL)"),
  A1C = list(file = "A1C.csv", ylabel = "Mean HbA1c (%)"),
  AST = list(file = "AST.csv", ylabel = "Mean AST (U/L)"),
  ALP = list(file = "ALP.csv", ylabel = "Mean ALP (U/L)")
)

# ------------------- Load MCB Patient Dataset -------------------
data_files <- c(
  "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv",
  "datasets/prob_validation_tapestry_with_predictions.csv",
  "datasets/tapestry_dbscan_assignments.csv",
  "datasets/centroid_mcb_test_centroid_assignments.csv"
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

# ------------------- Loop Over Measures -------------------
all_stats <- list()

for(measure in names(measure_config)) {
  config <- measure_config[[measure]]
  
  # Read measure file and merge with patient data
  data <- read.delim(config$file, sep = ",", header = TRUE) %>%
    merge(final_all_bio, ., by.x = "PATIENT_ID", by.y = "PATIENT_ID", all.x = TRUE)
  
  # Analyze measure
  results <- analyze_followup(
    data, 
    measure_col = "FLOWSHEET_RESULT_TXT", 
    file_prefix = measure,
    ylabel = config$ylabel
  )
  
  all_stats[[measure]] <- results$stat_test
}

# ------------------- Combine & Save All Statistics -------------------
combined_stats <- bind_rows(all_stats)
write.table(combined_stats, "all_measures_stats.csv", sep = "\t", row.names = FALSE, quote = FALSE)


setwd("follow_up/mcb/tapestry_fold")

measure_config <- list(
  BMI = list(file = "BMI.csv", ylabel = expression(paste("Mean BMI (kg/m"^{2},")"))),
  ALT = list(file = "ALT.csv", ylabel = "Mean ALT (U/L)"),
  HDL = list(file = "HDL.csv", ylabel = "Mean HDL (mg/dL)"),
  A1C = list(file = "A1C.csv", ylabel = "Mean HbA1c (%)"),
  AST = list(file = "AST.csv", ylabel = "Mean AST (U/L)"),
  ALP = list(file = "ALP.csv", ylabel = "Mean ALP (U/L)")
)

# ------------------- Load Tapestry Patient Dataset -------------------
data_files <- c(
  "datasets/prob_validation_tapestry_with_predictions.csv",
  "datasets/tapestry_dbscan_assignments.csv",
  "datasets/centroid_mcb_test_centroid_assignments.csv"
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

# ------------------- Loop Over Measures -------------------
all_stats <- list()

for(measure in names(measure_config)) {
  config <- measure_config[[measure]]
  
  # Read measure file and merge with patient data
  data <- read.delim(config$file, sep = ",", header = TRUE) %>%
    merge(final_all_bio, ., by.x = "PATIENT_ID", by.y = "PATIENT_ID", all.x = TRUE)
  
  # Analyze measure
  results <- analyze_followup(
    data, 
    measure_col = "FLOWSHEET_RESULT_TXT", 
    file_prefix = measure,
    ylabel = config$ylabel
  )
  
  all_stats[[measure]] <- results$stat_test
}

# ------------------- Combine & Save All Statistics -------------------
combined_stats <- bind_rows(all_stats)
write.table(combined_stats, "tap_all_measures_stats.csv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("✓ All analysis complete. Plots and statistics saved.\n")
