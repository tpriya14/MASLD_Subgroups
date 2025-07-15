# Install required packages if not already installed
pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")
all_pkgs <- c(pkgs_overall, pkgs_lca, "olsrr", "scatterplot3d", "MASS", "dplyr", "tidyr", "purrr", "stringr", "tidyfit", "table1")
lapply(all_pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    require(pkg, character.only = TRUE, quietly = TRUE)
  }
})

# Set working directory (update this path to your actual directory)
setwd("/MAYO_BIOBANK")

# Create results directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Read input data
final_all <- read.delim("l_info_02_06.txt", sep = "\t", header = TRUE)

# Select features for analysis
feature <- final_all %>% dplyr::select(
  AGE, PATIENT_GENDER, PATIENT_RACE_NAME, PATIENT_ETHNICITY_NAME, AVG_BMI,  # Sociodemographic and anthropometric
  ALT, AST, GGT, ALP, bilirubin,                                           # Liver-related lab tests
  A1C, GLUCOSE_FASTING, TRIGLYCERIDE, cholesterol, LDL, HDL, total_protein, # Other lab tests
  OBESITY, HYPERLIPIDEMIA, METABOLIC_SYNDROME, HEART_PROB_CODE, DIABETES    # Comorbidities
)

# Assess missing data
missing_data <- feature
missing_values <- is.na(missing_data)
percent_missing <- colMeans(missing_values) * 100
missing_df <- sort(percent_missing, decreasing = TRUE)
print(missing_df)

# Visualize missing data
par(mar = c(9, 9, 1, 1))
png("results/Missing_Value_Barplot.png", width = 1800, height = 1200, res = 200)
barplot(missing_df, 
        main = "Missingness Percentage of Clinical Lab Values",
        xlab = "Variables",
        ylab = "Missingness Percentage",
        col = "blue",
        ylim = c(0, 100),
        las = 3)  # Rotate x-axis labels vertically
dev.off()

# Convert comorbidity variables to Yes/No categories
final_all$Obesity1 <- ifelse(final_all$OBESITY == "", "No", "Yes")
final_all <- final_all %>% mutate(Hyperlipidemia1 = case_when(
  HYPERLIPIDEMIA == "" | is.na(HYPERLIPIDEMIA) ~ "No",
  TRUE ~ "Yes"
))
final_all <- final_all %>% mutate(MetS1 = case_when(
  METABOLIC_SYNDROME == "" | is.na(METABOLIC_SYNDROME) ~ "No",
  TRUE ~ "Yes"
))
final_all$CVD1 <- ifelse(final_all$HEART_PROB_CODE == "", "No", "Yes")
final_all$Sleep1 <- ifelse(final_all$SLEEP == "", "No", "Yes")
final_all$Hypertension1 <- ifelse(final_all$HYPERTENSION == "", "No", "Yes")
final_all$Joint1 <- ifelse(final_all$JOINT == "", "No", "Yes")
final_all$Depression1 <- ifelse(final_all$DEPRESSION == "", "No", "Yes")
final_all$Gerd1 <- ifelse(final_all$GERD == "", "No", "Yes")
final_all$Neurological1 <- ifelse(final_all$BRAIN == "", "No", "Yes")
final_all$Migraine1 <- ifelse(final_all$MIGRAINE == "", "No", "Yes")
final_all <- final_all %>% mutate(CKD1 = case_when(
  RENAL_DISEASE_PROBLEM == "" | is.na(RENAL_DISEASE_PROBLEM) ~ "No",
  TRUE ~ "Yes"
))


# Categorize BMI, Race, Ethnicity, and other clinical variables
final_all <- final_all %>% mutate(BMI_C = case_when(
  AVG_BMI >= 35 ~ "Obesity II/III",
  AVG_BMI >= 30 ~ "Obesity I",
  AVG_BMI >= 25 ~ "Overweight",
  AVG_BMI >= 18.5 ~ "Normal",
  AVG_BMI < 18.5 ~ "Underweight"
))

final_all <- final_all %>% mutate(PRS_C5 = case_when(
  PRS_5_SNP_SCORE1_SUM >= quantile(final_all$PRS_5_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
  PRS_5_SNP_SCORE1_SUM >= quantile(final_all$PRS_5_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(PRS_C10 = case_when(
  PRS_10_SNP_SCORE1_SUM >= quantile(final_all$PRS_10_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
  PRS_10_SNP_SCORE1_SUM >= quantile(final_all$PRS_10_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(PRS_C15 = case_when(
  PRS_15_SNP_SCORE1_SUM >= quantile(final_all$PRS_15_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
  PRS_15_SNP_SCORE1_SUM >= quantile(final_all$PRS_15_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(PRS_C68 = case_when(
  PRS_68_SNP_SCORE1_SUM >= quantile(final_all$PRS_68_SNP_SCORE1_SUM, probs = 0.9, na.rm = TRUE) ~ "High",
  PRS_68_SNP_SCORE1_SUM >= quantile(final_all$PRS_68_SNP_SCORE1_SUM, probs = 0.4, na.rm = TRUE) ~ "Intermediate",
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(GF_C = case_when(
  GLUCOSE_FASTING >= 126 ~ "Diabetes",
  GLUCOSE_FASTING >= 100 ~ "Prediabetes",
  is.na(GLUCOSE_FASTING) ~ NA_character_,
  TRUE ~ "Normal"
))
final_all <- final_all %>% mutate(A1C_C = case_when(
  A1C >= 6.5 ~ "Diabetes",
  A1C >= 5.7 ~ "Prediabetes",
  is.na(A1C) ~ NA_character_,
  TRUE ~ "Normal"
))
final_all <- final_all %>% mutate(ALT_C = case_when(
  ALT > 55 ~ "High",
  ALT >= 7 ~ "Normal",
  ALT < 7 ~ "Low",
  is.na(ALT) ~ NA_character_
))
final_all <- final_all %>% mutate(AST_C = case_when(
  AST > 48 ~ "High",
  AST >= 8 ~ "Normal",
  is.na(AST) ~ NA_character_,
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(ALP_C = case_when(
  ALP > 129 ~ "High",
  ALP >= 40 ~ "Normal",
  is.na(ALP) ~ NA_character_,
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(PT_C = case_when(
  PT > 12.5 ~ "High",
  PT >= 9.4 ~ "Normal",
  is.na(PT) ~ NA_character_,
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(BUN_C = case_when(
  BUN > 24 ~ "High",
  BUN >= 6 ~ "Normal",
  is.na(BUN) ~ NA_character_,
  TRUE ~ "Low"
))
final_all <- final_all %>% mutate(LDL_C = case_when(
  LDL < 0 | is.na(LDL) ~ NA_character_,
  LDL < 70 ~ "Best",
  LDL < 100 ~ "Optimal",
  LDL < 129 ~ "Optimal",
  LDL < 159 ~ "High",
  LDL < 189 ~ "High",
  TRUE ~ "High"
))
final_all <- final_all %>% mutate(HDL_C = case_when(
  HDL < 0 | is.na(HDL) ~ NA_character_,
  HDL < 50 ~ "Poor",
  HDL <= 59 ~ "Better",
  TRUE ~ "Best"
))
final_all <- final_all %>% mutate(TG_C =

 case_when(
  TRIGLYCERIDE < 0 | is.na(TRIGLYCERIDE) ~ NA_character_,
  TRIGLYCERIDE < 150 ~ "Desirable",
  TRIGLYCERIDE < 199 ~ "High",
  TRIGLYCERIDE < 499 ~ "High",
  TRIGLYCERIDE >= 500 ~ "High"
))
final_all <- final_all %>% mutate(Age_C = case_when(
  AGE < 45 ~ "< 45",
  AGE < 70 ~ "45 to 69",
  AGE >= 70 ~ "70 +"
))

# Summary of the dataset
summary(final_all)

# Save processed dataset
write.table(final_all, file = "T_MAYO_BIOBNAK_Pre_Process_Case_POOL_linical_info.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
