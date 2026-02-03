# =============================================
# Visualization Script for LCA Results
# =============================================
# This script generates figures comparing clinical variables (box plots)
# and disease prevalence (bar plots) across subgroups with statistical analysis 
# for three datasets:
#  - Development cohort (MCB 1A)
#  - Intenal validation cohort (MCB 1B)
#  - Independent validation cohort (Tapestry)
#
# Includes:
# - Boxplots with statistical comparisons for continuous variables
# - Grouped bar plots for disease prevalence
# - Statistical tests (Kruskal-Wallis, Dunn's test, Fisher's or Chi-square test)

# =============================================
# ------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------
required_packages <- c("dplyr","tidyr","purrr","ggplot2","rstatix","ggpubr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

run_visualization_pipeline <- function(dataset_name, input_file) {
  
  cat("\n========================================\n")
  cat("Processing Dataset:", dataset_name, "\n")
  cat("========================================\n\n")
  
  # ---- Create dataset-specific directories inside results/ ----
  dir.create(paste0("results/", dataset_name, "/figures/boxplots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("results/", dataset_name, "/figures/barplots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("results/", dataset_name, "/statistical_results"), showWarnings = FALSE, recursive = TRUE)
  
  
  cat("Loading data for visualization...\n")
  
  if(!file.exists(input_file)) {
    stop(paste("Error: Input file not found ->", input_file))
  }
  
  combined_data <- read.delim(
    input_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  cat("Data loaded:", nrow(combined_data), "samples\n")
  cat("Subgroup distribution:\n")
  print(table(combined_data$Subgroup))
  cat("\n")
  
  combined_data$Subgroup <- factor(
    combined_data$Subgroup,
    levels = c("Control", "C1", "C2", "C3", "C4", "C5")
  )
  
  
  # ------------------------------------------------------------
  # Function: Process Continuous Variables
  # ------------------------------------------------------------
  process_variable <- function(data, variable, subgroup = "Subgroup", adjust_method = "bonferroni") {
    
    cat(sprintf("Processing variable: %s\n", variable))
    
    data_clean <- data %>%
      filter(!is.na(.data[[variable]]), !is.na(.data[[subgroup]]))
    
    kruskal_result <- kruskal.test(
      as.formula(paste(variable, "~", subgroup)),
      data = data_clean
    )
    
    dunn_result <- data_clean %>%
      dunn_test(
        as.formula(paste(variable, "~", subgroup)),
        p.adjust.method = adjust_method
      ) %>%
      filter(group1 == "Control" | group2 == "Control")
    
    dunn_result <- dunn_result %>% add_xy_position(x = subgroup)
    dunn_result$kruskal_p <- kruskal_result$p.value
    
    # Set Y-axis label
    y_label <- switch(
      variable,
      "Age" = "Age (years)",
      "BMI" = "BMI (kg/m²)",
      "HDL" = "HDL (mg/dL)",
      "LDL" = "LDL (mg/dL)",
      "TRIGLYCERIDE" = "TG (mg/dL)",
      "PRS_15_SNP_SCORE1_SUM" = "PRS",
      "A1C" = "HbA1c (%)",
      paste(variable, "(U/L)")
    )
    
    # Y-axis limits
    y_limits <- switch(
      variable,
      "Age" = c(10, 80),
      "BMI" = c(10, 70),
      "HDL" = c(0, 100),
      "ALT" = c(0, 260),
      "AST" = c(0, 260),
      "ALP" = c(0, 260),
      "A1C" = c(2, 10),
      "PRS_15_SNP_SCORE1_SUM" = c(-2, 2),
      c(0, max(data_clean[[variable]], na.rm = TRUE) * 1.1)
    )
    
    step_increase_value <- ifelse(max(y_limits) > 100, 0.01, 0.05)
    tip_length <- ifelse(max(y_limits) > 100, 0.002, 0.01)
    
    # Create boxplot with previous colors, size, formatting
    p <- ggboxplot(
      data_clean, 
      x = subgroup, 
      y = variable, 
      fill = subgroup,
      bxp.errorbar = TRUE, 
      notch = TRUE,
      outlier.shape = 16,
      outlier.size = 1
    ) +
      scale_fill_manual(
        values = c(
          "#8B6914",      
          "#B0C4DE",      
          "#BBFFFF",     
          "darkseagreen1", 
          "lightsalmon",  
          "plum1"        
        )
      ) +
      coord_cartesian(ylim = y_limits * 1.9) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text = element_text(color = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none"
      ) +
      labs(y = y_label)
    
    # Set y-position for significance annotations
    y_position_value <- max(y_limits * 1.6)
    dunn_result <- dunn_result %>%
      mutate(y.position = y_position_value)
    
    # Add significance annotations
    p <- p + stat_pvalue_manual(
      dunn_result, 
      size = 6, 
      label = "p.adj.signif",
      tip.length = tip_length, 
      step.increase = step_increase_value,
      y.position = y_position_value
    )
    
    # Save plots in results/<dataset>/figures/boxplots
    ggsave(
      paste0("results/", dataset_name, "/figures/boxplots/", variable, ".png"), 
      plot = p, width = 6, height = 6, dpi = 400
    )
    
    ggsave(
      paste0("results/", dataset_name, "/figures/boxplots/", variable, ".pdf"), 
      plot = p, width = 6, height = 6
    )
    
    cat(sprintf("  Plot saved: results/%s/figures/boxplots/%s.png\n", dataset_name, variable))
    cat(sprintf("  Significant comparisons: %d\n\n", sum(dunn_result$p.adj < 0.05)))
    
    return(dunn_result)
  }
  
  # ------------------------------------------------------------
  # Analyze Continuous Variables
  # ------------------------------------------------------------
  cat("Analyzing Continuous Variables...\n")
  
  continuous_vars <- c("ALT", "AST", "ALP", "Age", "BMI", "A1C", "HDL")
  available_vars <- continuous_vars[continuous_vars %in% names(combined_data)]
  continuous_results <- list()
  
  for (variable in available_vars) {
    stat_results <- process_variable(combined_data, variable, "Subgroup")
    stat_results$Variable <- variable
    continuous_results[[variable]] <- stat_results
  }
  
  cat("Continuous analysis complete for:", dataset_name, "\n\n")
  
  
  # ------------------------------------------------------------
  # Disease Prevalence Analysis
  # ------------------------------------------------------------
  cat("Analyzing Disease Prevalence...\n")
  
  diseases <- c(
    "CKD1", "CVD1", "Depression1", "Diabetes1", "Hypertension1",
    "Hyperlipidemia1", "Joint1", "MetS1", "Migraine1",
    "Obesity1", "Sleep1"
  )
  available_diseases <- diseases[diseases %in% names(combined_data)]
  disease_data <- combined_data[, c(available_diseases, "Subgroup")]
  
  disease_data <- disease_data %>%
    mutate(across(all_of(available_diseases), as.factor))
  
  long_data <- disease_data %>%
    pivot_longer(
      cols = all_of(available_diseases),
      names_to = "Disease",
      values_to = "Status"
    ) %>%
    mutate(Status = ifelse(Status == "Yes" | Status == "1" | Status == 1, 1, 0))
  
  summary_data <- long_data %>%
    group_by(Subgroup, Disease) %>%
    summarise(
      Count = sum(Status),
      Total = n(),
      Proportion = Count / Total,
      SE = sqrt((Proportion * (1 - Proportion)) / Total),  # Standard error
      .groups = "drop"
    ) %>%
    mutate(
      Percentage = Proportion * 100,
      CI_lower = Percentage - (1.96 * SE * 100),  # 95% CI
      CI_upper = Percentage + (1.96 * SE * 100)
    )
  
  # Barplot
  p_disease <- ggplot(summary_data, aes(x = Disease, y = Percentage, fill = Subgroup)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, color = "gray40") +
    geom_errorbar( aes(ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.9), width = 0.2, color = "black", size = 0.3 ) +
    scale_fill_manual(
      values = c(
        "#8B6914", "#B0C4DE", "#BBFFFF", "darkseagreen1", "lightsalmon", "plum1"
      )
    ) +
    labs(x = "", y = "Prevalence (%)", fill = "Subgroup") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
      panel.background = element_rect(fill = "white"),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Save barplots in results/<dataset>/figures/barplots
  ggsave(
    paste0("results/", dataset_name, "/figures/barplots/disease_prevalence.png"),
    plot = p_disease, width = 12, height = 7, dpi = 400
  )
  
  ggsave(
    paste0("results/", dataset_name, "/figures/barplots/disease_prevalence.pdf"),
    plot = p_disease, width = 12, height = 7
  )
  
  write.csv(
    summary_data,
    paste0("results/", dataset_name, "/statistical_results/statistical_disease_prevalence_summary.csv"),
    row.names = FALSE
  )
  
  cat("Disease prevalence analysis complete for:", dataset_name, "\n\n")
}


# ------------------------------------------------------------
# RUN PIPELINE FOR ALL THREE DATASETS
# ------------------------------------------------------------
run_visualization_pipeline(
  "development",
  "datasets/development_data_for_visualization.csv"
)

run_visualization_pipeline(
  "test",
  "datasets/test_data_for_visualization.csv"
)

run_visualization_pipeline(
  "tapestry",
  "datasets/tapestry_test_data_for_visualization.csv"
)

cat("\n========================================\n")
cat("ALL VISUALIZATIONS COMPLETE!\n")
cat("========================================\n")
