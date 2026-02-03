# =============================================
# MRE Analysis
# =============================================
# This script evaluates the distribution of MRE (magnetic resonance elastography) measurements 
# across clinically defined subgroups in both development and validation cohorts.
#
# Objectives:
# - Assess MRE proportions within subgroups in the MCB development cohort and its validation cohort
# - Evaluate MRE distributions in the independent Tapestry dataset using centroid-based subgroup assignment methods
# =============================================

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

process_mre_dataset <- function(input_file, output_prefix, id_var) {
  data <- read.delim(input_file, sep = "\t", header = TRUE)
  MRE1 <- data
  
  # ------------------- Categorize Stiffness -------------------
  MRE1 <- MRE1 %>% mutate(
    fibro_pbm = case_when(
      impression_text < 2.5 ~ 'LS<2.5',
      impression_text <= 5 ~ 'LS(2.5 to 5)',
      impression_text > 5 ~ 'LS>5',
      TRUE ~ as.character(impression_text)
    )
  )

  merged_data <- merge(data, MRE1[, c(1,5,10,11)], 
                       by.x = id_var, by.y = id_var, all.x = TRUE)
  
  percentage_data <- merged_data %>%
    group_by(fibro_pbm, Subgroup) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Subgroup) %>%
    mutate(
      total_count = sum(Count),
      Percentage = (Count / sum(Count)) * 100
    ) %>%
    filter(!is.na(fibro_pbm))
  
  percentage_data$fibro_pbm <- factor(
    percentage_data$fibro_pbm, 
    levels = c("LS<2.5", "LS(2.5 to 5)", "LS>5")
  )
  
  # ------------------- Create Plot -------------------
  p <- ggplot(percentage_data, aes(x = fibro_pbm, y = Percentage, fill = Subgroup)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1), width = 0.8) +
    geom_text(aes(label = paste0(Count, '(', round(Percentage, 1), '%)')),
              position = position_dodge(width = 1),
              vjust = 0.5, hjust = -0.1, size = 4.5, angle = 90) +
    labs(y = "Percentages (%)") +
    ylim(0,50) +
    scale_fill_manual(values = c("#B0C4DE", "#BBFFFF", "darkseagreen1", 
                                 "lightsalmon", "plum1")) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 20, color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 16, color = "black")
    )

  print(p)
  
  ggsave(paste0(output_prefix, "_MRE_plot.png"), plot = p, width = 7, height = 5, dpi = 500)
  ggsave(paste0(output_prefix, "_MRE_plot.pdf"), plot = p, width = 7, height = 5, dpi = 500)
  write.csv(percentage_data, paste0(output_prefix, "_MRE_summary.csv"), row.names = FALSE)
  
  return(percentage_data)
}

# ------------------- Run Analysis for Development Cohort -------------------
bio_results <- process_mre_dataset(
  input_file = "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv",
  output_prefix = "bio",
  id_var = "PATIENT_ID"
)

# ------------------- Run Analysis for Validation Cohort -------------------
tap_results <- process_mre_dataset(
  input_file = "datasets/centroid_mcb_test_centroid_assignments.csv",
  output_prefix = "tap",
  id_var = "PATIENT_ID"
)

# ------------------- Combine Results -------------------
combined <- rbind(bio_results, tap_results)

combined_summary <- combined %>%
  group_by(fibro_pbm, Subgroup) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Subgroup) %>%
  mutate(
    Percentage = (Count / sum(Count)) * 100
  )

write.csv(combined_summary, "combined_MRE_summary.csv", row.names = FALSE)

