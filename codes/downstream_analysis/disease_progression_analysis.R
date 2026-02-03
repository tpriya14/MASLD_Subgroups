# =============================================
# Disease Progression Analysis
# =============================================
# This script analyzes intrahepatic and extrahepatic disease 
# progression risk across subgroups using Cox proportional 
# hazards models and cumulative incidence plots.
#
# Objectives:
# - Evaluate Cox regression and cumulative incidence in MCB 
#   development and validation cohorts
# - Evaluate Cox regression and cumulative incidence in the 
#   Tapestry dataset using three assignment methods
# =============================================

# ---------------------- Load Libraries ---------------------- #
required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr",
  "ggplot2", "survival", "survminer", "tibble", "lubridate",
  "ggsurvfit", "tidycmprsk", "gtsummary", "readr",
  "caret", "randomForest", "Boruta", "mlbench",
  "reshape2", "ggrepel", "scales", "paletteer",
  "poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

invisible(lapply(required_packages, library, character.only = TRUE))

# Function to prepare data, run Cox model, summarize, and plot survival
analyze_survival <- function(data, PATIENT_ID_col = "PATIENT_ID", status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix) {
  # Filter and clean data
  
  # Year1 represents the time difference (in years) between the MASLD diagnosis date
  # and the first occurrence of the event of interest.
  # The follow-up period is up to 10 years after MASLD diagnosis.
  # Any events occurring before MASLD diagnosis or at time 0 are excluded.
  
  max_year = 10
  filtered_data <- data %>%
    filter(!is.na(Date)) %>%
    group_by(across(all_of(c(PATIENT_ID_col, subgroup_col)))) %>%
    filter(!(Year1 <= 0 & !!sym(status_col) == 1)) %>%
    filter(Year1 >= 0 & Year1 <= max_year)
  
  # Cox proportional hazards model
  fit <- coxph(Surv(Year1, !!sym(status_col)) ~ !!sym(subgroup_col), data = filtered_data)

  summary(fit)
  cox_summary <- fit
  # Function to fit Cox model and extract summary
  cox_results <- tidy(cox_summary) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = gsub("Subgroup", "", term) ,
      HR = exp(estimate),
      lower_95 = exp(estimate - 1.96 * std.error),
      upper_95 = exp(estimate + 1.96 * std.error),
      p.value= paste0(sprintf("%.4f", p.value)),
      HR_CI = paste0(sprintf("%.2f", HR), "(", sprintf("%.2f", lower_95), "-", sprintf("%.2f", upper_95), ")")
    ) %>%
    select(term, HR_CI, p.value)
  cox_results
  default_value <- data.frame(
    term= "C1",
    HR_CI = "1.00",  # Example HR_CI
    p.value = '-'  # Example p.value
  )
  default_value
  summary_table <- rbind( default_value,cox_results)
  summary_table
  names(summary_table) <- c("Subgroup", "HR_CI", "p.value")
  
  # Calculate the number of events and censored cases for each subgroup
  
  filtered_data$Status <- as.factor(filtered_data$Status)
  event_data <- filtered_data %>%
    group_by(subgroup_col) %>%
    summarize(
      Events = sum(Status == 1),
      Censored = sum(Status == 0),
      PctEvents = Events / (Events + Censored) * 100
    )
  #View(date_diff_data)
  event_data
  cox_results
  
  # Combine the Cox results with event data
  summary_table1 <-summary_table %>%
    left_join(event_data, by = "Subgroup") %>%
    select(subgroup_col, Events, Censored, PctEvents, HR_CI, p.value)
  
  # Rename columns for better readability
  summary_table <- summary_table1 %>%
    rename(
      `Subgroup` = subgroup_col,
      `Events` = Events,
      `Censored` = Censored,
      `% Events` = PctEvents,
      `HR (95% CI)` = HR_CI,
      `P Value` = p.value
    )
  summary_table
  
  # Save summary table
  write.table(summary_table, paste0(plot_file_prefix, "_survival_summary.csv"), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Survival plot
  surv_plot <- cuminc(Surv(Year1, !!sym(status_col)) ~ !!sym(subgroup_col), data = filtered_data) %>%
    ggcuminc(linewidth=1) +
    labs(
      x = "Years"
    ) +
    scale_x_continuous(breaks=seq(from = 0, to =10, by = 2))+
    theme(
      panel.background = element_rect(fill = "white"),  # set background to white
      panel.grid.major.x = element_blank(),  # remove major gridlines on y-axis
      panel.grid.major.y = element_blank(),  # remove major gridlines on y-axis
      panel.grid.minor.x = element_blank(),  # remove major gridlines on y-axis
      panel.grid.minor.y = element_blank(),  # remove minor gridlines on y-axis
      axis.line.x = element_line(color = "black"),  # keep x-axis line
      axis.line.y = element_line(color = "black"),  # keep y-axis line
      axis.line = element_line(size = 0.5),  # customize axis line size
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "none"
    ) +
    add_risktable(
      risktable_stats = "n.risk",
      stats_label = list(n.risk = "Number at Risk"),
      size = 4.8,  # Increase the size of the risk table text
      risktable_height = .25,
      theme = theme_risktable_default(axis.text.y.size = 13, plot.title.size = 14)
    ) 
  surv_plot <- surv_plot + 
    scale_color_manual(values = c("steelblue", "turquoise", "seagreen", "salmon", "purple"))
  
  surv_plot
  # Save plots
  ggsave(paste0(plot_file_prefix, "_survival_plot.png"), plot = surv_plot, width = 5, height = 6, dpi = 500)
  ggsave(paste0(plot_file_prefix, "_survival_plot.pdf"), plot = surv_plot, width = 5, height = 6, dpi = 500)
  
  return(list(summary_table = summary_table, surv_plot = surv_plot))
}

# Load the data files
ckd_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_ckd_assignments.csv", sep = "\t", header = TRUE)
cvd_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_cvdd_assignments.csv", sep = "\t", header = TRUE)
sleep_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_sleep_assignments.csv", sep = "\t", header = TRUE)
depression_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_depression_assignments.csv", sep = "\t", header = TRUE)
mash_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_mash_assignments.csv", sep = "\t", header = TRUE)
fibrosis_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_fibrosis_assignments.csv", sep = "\t", header = TRUE)
cirrhosis_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_cirrhosis_assignments.csv", sep = "\t", header = TRUE)
hcc_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_hcc_assignments.csv", sep = "\t", header = TRUE)
lt_date_diff_data1 <- read.delim("follow_up/mcb/followup_mcb_centroid_lt_assignments.csv", sep = "\t", header = TRUE)


ckd_results <- analyze_survival(ckd_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "ckd_failure")
cvd_results <- analyze_survival(cvd_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "CVD_ischemic")
sleep_results <- analyze_survival(sleep_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "sleep")
depression_results <- analyze_survival(depression_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "depressive_major")
mash_results <- analyze_survival(mash_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "MASH")
fibrosis_results <- analyze_survival(fibrosis_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "Fibrosis")
cirrhosis_results <- analyze_survival(cirrhosis_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "Cirrhosis")
hcc_results <- analyze_survival(hcc_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "HCC")
lt_results <- analyze_survival(lt_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "LT")

# Repeat the same analysis with the Tapestry dataset with different assignmnet methods
# Load the data files
ckd_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_ckd_assignments.csv", sep = "\t", header = TRUE)
cvd_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_cvdd_assignments.csv", sep = "\t", header = TRUE)
sleep_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_sleep_assignments.csv", sep = "\t", header = TRUE)
depression_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_depression_assignments.csv", sep = "\t", header = TRUE)
mash_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_mash_assignments.csv", sep = "\t", header = TRUE)
fibrosis_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_fibrosis_assignments.csv", sep = "\t", header = TRUE)
cirrhosis_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_cirrhosis_assignments.csv", sep = "\t", header = TRUE)
hcc_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_hcc_assignments.csv", sep = "\t", header = TRUE)
lt_date_diff_data1 <- read.delim("follow_up/tapestry/followup_tapestry_centroid_lt_assignments.csv", sep = "\t", header = TRUE)


ckd_results <- analyze_survival(ckd_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_ckd_failure")
cvd_results <- analyze_survival(cvd_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_CVD_ischemic")
sleep_results <- analyze_survival(sleep_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_sleep")
depression_results <- analyze_survival(depression_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_depressive_major")
mash_results <- analyze_survival(mash_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_MASH")
fibrosis_results <- analyze_survival(fibrosis_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_Fibrosis")
cirrhosis_results <- analyze_survival(cirrhosis_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_Cirrhosis")
hcc_results <- analyze_survival(hcc_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_HCC")
lt_results <- analyze_survival(lt_date_diff_data1, status_col = "Status", subgroup_col = "Subgroup", plot_file_prefix = "tap_LT")
