# =============================================
# Data Preprocessing Script
# =============================================
# This script loads MCB and Tapestry raw data and creates categorical
# variables according to different cutoff values.
# =============================================

# ------------------- Load Required Packages -------------------
required_packages <- c("dplyr", "tidyr", "table1")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------- Define Dataset Processing Function -------------------
process_MASLD_dataset <- function(
    input_file,
    output_file,
    dataset_name = "Dataset",
    make_table1 = FALSE,
    table1_data = NULL,
    table1_output = NULL
) {
  
  # ------------------- Load Data -------------------
  final_all <- read.delim(
    input_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # ------------------- BMI -------------------
  final_all <- final_all %>%
    mutate(
      BMI_C = case_when(
        BMI >= 35 ~ "Obesity II/III",
        BMI >= 30 ~ "Obesity I",
        BMI >= 25 ~ "Overweight",
        BMI >= 18.5 ~ "Normal",
        BMI < 18.5 ~ "Underweight",
        TRUE ~ NA_character_
      )
    )
  final_all$BMI_C <- factor(
    final_all$BMI_C,
    levels = c("Underweight", "Normal", "Overweight", "Obesity I", "Obesity II/III")
  )
  
  # ------------------- PRS Categories -------------------
  prs_cut <- function(x) {
    case_when(
      x >= quantile(x, 0.9, na.rm = TRUE) ~ "High",
      x >= quantile(x, 0.4, na.rm = TRUE) ~ "Intermediate",
      TRUE ~ "Low"
    )
  }
  final_all <- final_all %>%
    mutate(
      PRS_C15 = prs_cut(PRS_15_SNP_SCORE1_SUM),
      PRS_C5  = prs_cut(PRS_5_SNP_SCORE1_SUM)
    )
  prs_levels <- c("Low", "Intermediate", "High")
  final_all$PRS_C15 <- factor(final_all$PRS_C15, levels = prs_levels)
  final_all$PRS_C5  <- factor(final_all$PRS_C5,  levels = prs_levels)
  
  # ------------------- Glucose / A1C -------------------
  final_all <- final_all %>%
    mutate(
      A1C_C = case_when(
        A1C >= 6.5 ~ "Diabetes",
        A1C >= 5.7 ~ "Prediabetes",
        TRUE ~ "Normal"
      )
    )
  final_all$A1C_C <- factor(final_all$A1C_C, levels = c("Normal", "Prediabetes", "Diabetes"))
  
  # ------------------- Liver Enzymes -------------------
  enzyme_levels <- c("Low", "Normal", "High")
  final_all <- final_all %>%
    mutate(
      ALT_C = case_when(
        ALT > 55 ~ "High",
        ALT >= 7 ~ "Normal",
        TRUE ~ "Low"
      ),
      AST_C = case_when(
        AST > 48 ~ "High",
        AST >= 8 ~ "Normal",
        TRUE ~ "Low"
      ),
      ALP_C = case_when(
        ALP > 129 ~ "High",
        ALP >= 40 ~ "Normal",
        TRUE ~ "Low"
      )
    )
  final_all$ALT_C <- factor(final_all$ALT_C, levels = enzyme_levels)
  final_all$AST_C <- factor(final_all$AST_C, levels = enzyme_levels)
  final_all$ALP_C <- factor(final_all$ALP_C, levels = enzyme_levels)
  
  # ------------------- BUN -------------------
  final_all <- final_all %>%
    mutate(
      BUN_C = case_when(
        BUN > 24 ~ "High",
        BUN >= 6 ~ "Normal",
        TRUE ~ "Low"
      )
    )
  final_all$BUN_C <- factor(final_all$BUN_C, levels = enzyme_levels)
  
  # ------------------- Lipids -------------------
  final_all <- final_all %>%
    mutate(
      LDL_C = case_when(
        LDL < 70 ~ "Best",
        LDL < 129 ~ "Optimal",
        TRUE ~ "High"
      ),
      HDL_C = case_when(
        HDL < 50 ~ "Poor",
        HDL <= 59 ~ "Better",
        TRUE ~ "Best"
      ),
      TG_C = case_when(
        TRIGLYCERIDE < 150 ~ "Desirable",
        TRUE ~ "High"
      )
    )
  final_all$LDL_C <- factor(final_all$LDL_C, levels = c("Best", "Optimal", "High"))
  final_all$HDL_C <- factor(final_all$HDL_C, levels = c("Poor", "Better", "Best"))
  final_all$TG_C  <- factor(final_all$TG_C,  levels = c("Desirable", "High"))
  
  # ------------------- Age -------------------
  final_all <- final_all %>%
    mutate(
      Age_C = case_when(
        Age < 45 ~ "< 45",
        Age < 70 ~ "45 to 69",
        TRUE ~ "70 +"
      )
    )
  final_all$Age_C <- factor(final_all$Age_C, levels = c("< 45", "45 to 69", "70 +"))
  
  # ------------------- Save Processed Data -------------------
  write.table(final_all, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # ------------------- Table 1 -------------------
  if (make_table1 && !is.null(table1_data) && !is.null(table1_output)) {
    my.render.cat <- function(x) {
      c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
    }
    exclude_cols <- c("PATIENT_ID")
    table_columns <- table1_data[, !(names(table1_data) %in% exclude_cols)]
    
    x <- table1(
      ~ . | MASLD,
      data = table_columns,
      overall = TRUE,
      render.categorical = my.render.cat,
      topclass = "Rtable1-zebra"
    )
    
    tab1_df <- as.data.frame(x)
    
    write.table(tab1_df, table1_output, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  return(final_all)
}

# ------------------- Process MCB Dataset -------------------
MCB_processed <- process_MASLD_dataset(
  input_file  = "/MCB_Data.txt",
  output_file = "/MCB_Data_Processed.txt",
  dataset_name = "MCB",
  make_table1 = TRUE,
  table1_data = LCA_MASLD_data_75,
  table1_output = "/datasets/summary_by_MASLD_MCB.csv"
)

# ------------------- Process Tapestry Dataset -------------------
Tapestry_processed <- process_MASLD_dataset(
  input_file  = "/Tapestry_Data.txt",
  output_file = "/Tapestry_Data_Processed.txt",
  dataset_name = "Tapestry",
  make_table1 = TRUE,
  table1_data = LCA_MASLD_data_75,
  table1_output = "/datasets/summary_by_MASLD_Tapestry.csv"
)
