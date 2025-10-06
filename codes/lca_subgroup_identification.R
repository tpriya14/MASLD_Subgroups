# Latent Class Analysis for MASLD Biobank Data (Cohort 1A)
# =============================================
# This script performs Latent Class Analysis (LCA) on Metabolic dysfunction–Associated Steatotic Liver Disease (MASLD) biobank data
# to identify latent clusters based on demographic and clinical variables. The analysis is conducted on both
# the full dataset (100%) and a 75% training subset, with model performance evaluated using AIC, BIC, SABIC,
# entropy, and likelihood ratio tests. Results are saved as publication-ready tables and datasets.


# Load Required Packages
# ----------------------
library(data.table)       # Efficient data handling
library(dplyr)            # Data manipulation
library(tidyr)            # Data tidying
library(purrr)            # Functional programming
library(ggplot2)          # Visualization
library(poLCA)            # Latent Class Analysis
library(MASS)             # Statistical functions
library(olsrr)            # Stepwise regression
library(table1)           # Summary tables
library(tidyLPA)          # Tidy Latent Profile Analysis
library(rstatix)          # For Dunn's test and statistical annotations
library(ggpubr)           # For creating boxplots with statistical annotations

# Set Working Directory
# ---------------------
# Note: Update this path to your local directory or use relative paths for portability.
setwd("/MAYO_BIOBANK/")

# Load Data
# ---------
# Load pre-processed dataset with baseline demographic and clinical variables.
final_all <- read.delim("MCB_Data.txt", 
                        sep = "\t", header = TRUE)

# Data Preparation
# ----------------
# Remove unnecessary columns and filter for MASLD cases.
feature <- final_all
LCA_MASLD_data <- feature[feature$MASLD == 1, ]
selected_data_100 <- na.omit(LCA_MASLD_data)  # Remove rows with missing values

# Select variables for LCA and convert character columns to factors.
LCA_data <- selected_data_100
character_columns <- sapply(LCA_data, is.character)
LCA_data[, character_columns] <- lapply(LCA_data[, character_columns], as.factor)

# Set Seed for Reproducibility
# ----------------------------
set.seed(5055)  # Ensures consistent results across runs

# Split Data into Training (75%) and Test (25%) Sets
# --------------------------------------------------
test_index <- sample(1:nrow(selected_data_100), nrow(selected_data_100) * 0.25)
LCA_MASLD_data_75 <- selected_data_100[-test_index, ]

# Prepare Control Data (MASLD == 0) for Feature Selection
# -------------------------------------------------------
LCA_CTRL_data <- feature[feature$MASLD == 0, ]
selected_data_ctrl <- LCA_CTRL_data[, -c(1:5, 7:29, 44, 48, 49, 51, 61:65)]
selected_data_ctrl <- na.omit(selected_data_ctrl)

# Combine Training and Control Data for Feature Selection
# -------------------------------------------------------
feature_all1 <- rbind(LCA_MASLD_data_75[, names(LCA_data)], selected_data_ctrl)
feature_all1[, character_columns] <- lapply(feature_all1[, character_columns], as.numeric)

# Feature Selection Using Stepwise Regression
# -------------------------------------------
# Note: Stepwise regression is used here, but consider justifying or exploring alternatives (e.g., LASSO)
# for robustness in the manuscript.
model <- glm(MASLD ~ ., data = feature_all1)
step_both <- stepAIC(model, direction = "both")
summary(step_both)

# Define Variables for LCA
# ------------------------
# Selected based on stepwise regression results
vari_lca <- c("PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", "MetS1", "Hypertension1", 
              "ALT_C", "AST_C", "BMI_C", "HDL_C", "Depression1", "Migraine1", "CKD1", "ALP_C")

# Prepare LCA Data for 75% Dataset
# ----------------------------------------
dat_lca_100 <- LCA_MASLD_data_75[, vari_lca]

# Latent Class Analysis (LCA) Model Fitting
# -----------------------------------------
# Define LCA formula
f_lca <- as.formula(paste0("cbind(", paste(vari_lca, collapse = ", "), ") ~ 1"))

# Fit LCA models with 1 to 10 classes
n_start <- 1
n_end <- 10
out_path <- "/result/"
for (i in n_start:n_end) {
  mi_start <- poLCA(f_lca, data = dat_lca_100, nclass = i, maxiter = 1000, nrep = 30, 
                    graph = FALSE, calc.se = TRUE, verbose = FALSE)
  vprobs.start <- poLCA.reorder(mi_start$probs.start, order(mi_start$P, decreasing = TRUE))
  mi <- poLCA(f_lca, data = dat_lca_100, nclass = i, maxiter = 1000, 
              probs.start = vprobs.start, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  save(mi_start, mi, file = paste0(out_path, "lca_m", i, ".RData"))
}

# Save LCA data for reproducibility
save(dat_lca_100, file = paste0(out_path, "dat_lca_100.RData"))

# Evaluate LCA Model Performance
# ------------------------------
res_lca <- data.table()
for (j in n_start:n_end) {
  load(paste0(out_path, "lca_m", j, ".RData"))
  prop <- table(mi$predclass) / mi$Nobs * 100
  pre_prop <- mi$P * 100
  val_post <- apply(mi$posterior, 1, max)
  sum_post <- summary(val_post)
  
  mean_post <- if (j == 1) 1 else {
    dat_post <- data.table(posterior = val_post, class = mi$predclass)
    dat_post[, mean(posterior), by = "class"]$V1
  }
  
  nadj <- (mi$Nobs + 2) / 24
  sabic <- -2 * mi$llik + mi$npar * log(nadj)
  
  entropy <- function(p) sum(-p * log(p))
  error_prior <- entropy(mi$P)
  error_post <- mean(apply(mi$posterior, 1, entropy), na.rm = TRUE)
  entropy_val <- round(((error_prior - error_post) / error_prior), 3)
  
  res_j <- data.table(
    n.class = j, sample.size = mi$Nobs, n.param = mi$npar, loglike = mi$llik,
    AIC = mi$aic, BIC = mi$bic, SABIC = sabic, entropy = entropy_val,
    min.post = sum_post[1], median.post = sum_post[3], mean.post = sum_post[4],
    max.post = sum_post[6], mean.per.post = paste(round(mean_post, 2), collapse = "|"),
    min.per.post = min(mean_post), pre.prop = paste(round(pre_prop, 2), collapse = "|"),
    min.preprop = min(pre_prop), class.prop = paste(round(prop, 2), collapse = "|"),
    min.prop = min(prop)
  )
  res_lca <- rbind(res_lca, res_j)
}

# Calculate p-values for Likelihood Ratio Tests (LRT)
for (n in 2:n_end) {
  plrt <- calc_lrt(n = res_lca$sample.size[n], null_ll = res_lca$loglike[n-1], 
                   null_param = res_lca$n.param[n-1], null_classes = res_lca$n.class[n-1], 
                   alt_ll = res_lca$loglike[n], alt_param = res_lca$n.param[n], 
                   alt_classes = res_lca$n.class[n])[4]
  res_lca$p_lrt[n] <- plrt
}

# Save model performance results
write.table(res_lca, "res_lca_100.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# Fit Final LCA Model with 5 Classes
# ----------------------------------
mi_start <- poLCA(f_lca, data = dat_lca_100, nclass = 5, maxiter = 1000, nrep = 30, 
                  graph = FALSE, calc.se = TRUE, verbose = FALSE)
vprobs.start <- poLCA.reorder(mi_start$probs.start, order(mi_start$P, decreasing = TRUE))
lc51 <- poLCA(f_lca, data = dat_lca_100, nclass = 5, maxiter = 1000, 
              probs.start = vprobs.start, graph = FALSE, calc.se = TRUE, verbose = FALSE)

# Assign latent class clusters
LCA_MASLD_data_75$LatentClassCluster <- lc51$predclass

# Define custom rendering for categorical variables
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

# Create publication-ready summary table
table_column <- LCA_MASLD_data_75[, -c(1:5, 7, 21, 23, 25, 27, 29)]
x <- table1(~ . | LatentClassCluster, data = table_column, overall = TRUE, 
            render.categorical = my.render.cat, topclass = "Rtable1-zebra")

# Save summary table and clustered data
tab1_df <- as.data.frame(x)
write.table(tab1_df, "Bio_After_LCA_100.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(LCA_MASLD_data_75, "T100_5_Tapestry_After_LCA_T2043_clinical_info_06_11.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Follow the same steps as for the whole dataset. Fit and evaluate LCA models for 100% dataset to use the resutls as reference result

# Perform statistical analysis and visualization of clinical variables and disease prevalence
#          across subgroups in the MAYO Biobank dataset (Cohort 1A).

selected_data_75 <- LCA_MASLD_data_75
# Levels ensure consistent ordering: Control, C1, C2, C3, C4, C5
selected_data_75$Subgroup <- factor(selected_data_75$Subgroup, 
                                   levels = c("Control", "C1", "C2", "C3", "C4", "C5"))


# Print column names to verify data structure
colnames(selected_data_75)

# Function: process_variable
# Purpose: Perform Kruskal-Wallis and Dunn's tests for a given variable across subgroups,
#          and generate boxplots with significance annotations
# Arguments:
#   - data: Input data frame
#   - variable: Name of the variable to analyze (e.g., "ALT", "Age")
#   - subgroup: Name of the subgroup column (default: "Subgroup")
#   - adjust_method: P-value adjustment method (default: "bonferroni")
# Returns: Data frame with Dunn's test results
process_variable <- function(data, variable, subgroup, adjust_method = "bonferroni") {
  # Ensure subgroup is a factor with correct levels
  data[[subgroup]] <- factor(data[[subgroup]], levels = c("Control", "C1", "C2", "C3", "C4", "C5"))
  
  # Print subgroup column name for debugging
  print(subgroup)
  
  # Perform Kruskal-Wallis test to compare variable across subgroups
  kruskal_test_result <- kruskal.test(as.formula(paste(variable, "~", subgroup)), data = data)
  print(kruskal_test_result)
  
  # Perform Dunn's test for pairwise comparisons, focusing on Control vs. others
  dunn_result <- data %>%
    dunn_test(as.formula(paste(variable, "~", subgroup)), p.adjust.method = adjust_method) %>%
    filter(group1 == "Control" | group2 == "Control")
  
  # Add x-y positions for significance annotations in plots
  dunn_result <- dunn_result %>% add_xy_position(x = "Subgroup")
  
  # Add Kruskal-Wallis p-value to Dunn's test results
  dunn_result$kruskal_p <- kruskal_test_result$p.value
  
  # Format adjusted p-values for display
  dunn_result$p.format <- p_format(
    dunn_result$p.adj, accuracy = 0.05,
    leading.zero = FALSE
  )
  
  # Define y-axis label based on variable
  y_label <- switch(variable,
                    "Age" = "Age (years)",
                    "BMI" = "BMI (kg/m²)",
                    "HDL" = "HDL (mg/dL)",
                    "LDL" = "LDL (mg/dL)",
                    "TRIGLYCERIDE" = "TG (mg/dL)",
                    "PRS_15_SNP_SCORE1_SUM" = "PRS",
                    "A1C" = "A1C (%)",
                    paste(variable, "(U/L)"))  # Default for other variables
  
  # Define y-axis limits based on variable
  y_limits <- switch(variable,
                     "Age" = c(10, 80),
                     "BMI" = c(10, 70),
                     "HDL" = c(0, 100),
                     "ALT" = c(0, 260),
                     "AST" = c(0, 260),
                     "ALP" = c(0, 260),
                     "A1C" = c(2, 10),
                     c(0, max(data[[variable]], na.rm = TRUE) * 1.1))  # Default based on max value
  
  # Print Dunn's test results for verification
  print(dunn_result)
  
  # Set step increase and tip length for significance annotations
  step_increase_value <- ifelse(max(y_limits) > 100, 0.01, 0.05)
  tip_length <- ifelse(max(y_limits) > 100, 0.002, 0.01)
  
  # Create boxplot with custom styling
  p <- ggboxplot(data, x = subgroup, y = variable, fill = subgroup, 
                 bxp.errorbar = TRUE, notch = TRUE) +
    scale_fill_manual(values = c("#8B6914", "#B0C4DE", "#BBFFFF", 
                                 "darkseagreen1", "lightsalmon", "plum1")) +
    coord_cartesian(ylim = y_limits * 1.9) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 15),
      axis.line = element_line(size = 0.5)
    ) +
    theme(legend.position = "none") +
    labs(x = "Subgroups", y = y_label) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 15)
    )
  
  # Print Dunn's test results again for verification
  print(dunn_result)
  
  # Dynamically set y-position for significance annotations (90% of max y-limit)
  y_position_value <- max(y_limits * 1.6)
  dunn_result <- dunn_result %>%
    mutate(y.position = y_position_value)
  
  # Add significance annotations to the plot
  p <- p + stat_pvalue_manual(dunn_result, size = 6, label = "p.adj.signif",
                              tip.length = tip_length, step.increase = step_increase_value, 
                              y.position = y_position_value)
  
  # Save plot as PNG and PDF
  ggsave(paste(variable, ".png", sep = ""), plot = p, width = 4, height = 6, dpi = 400)
  ggsave(paste(variable, ".pdf", sep = ""), plot = p, width = 4, height = 6, dpi = 400)
  
  # Print Dunn's test results for final verification
  print(dunn_result)
  
  return(dunn_result)
}

# Define list of continuous variables to analyze
variables <- c("ALT", "AST", "ALP", "Age", "BMI", "A1C", "HDL")

# Initialize empty tibble to store statistical results
combined_results <- tibble()

# Loop through variables and process each one
for (variable in variables) {
  # Perform statistical tests and generate boxplots
  stat.test <- process_variable(selected_data_75, variable, "Subgroup")
  
  # Add variable name to results
  stat.test <- stat.test %>%
    mutate(Variable = variable)
  
  # Combine results into a single tibble
  combined_results <- bind_rows(combined_results, stat.test)
}

# Convert list columns to strings for CSV export
combined_results <- combined_results %>%
  mutate(across(where(is.list), ~sapply(., toString)))

# Export combined statistical results to CSV
write.csv(combined_results, "p_lab_combined_bio_clinical_stat_results.csv", row.names = FALSE)

# Print column names of the dataset for verification
colnames(selected_data_75)

# Define list of disease variables to analyze
diseases <- c("CKD1", "CVD1", "Depression1", "Diabetes1", "Hypertension1",
              "Hyperlipidemia1", "Joint1", "MetS1", "Migraine1",
              "Obesity1", "Sleep1", "Transplant", "Fibrosis", 
              "Cirrhosis", "Hepatocellular", "Steatohepatitis")

# Subset data to include only disease variables and Subgroup
data_file_1 <- selected_data_75[, c(diseases, "Subgroup")]

# Convert Diabetes1 to binary (Yes/No) for consistency
data_file_1 <- data_file_1 %>%
  mutate(Diabetes1 = ifelse(Diabetes1 == "Type 2", "Yes", "No"))

# Update selected_data_75 with processed data
selected_data_75 <- data_file_1

# Convert Steatohepatitis to factor
selected_data_75$Steatohepatitis <- as.factor(selected_data_75$Steatohepatitis)

# Print summary of the dataset
summary(selected_data_75)

# Reshape data to long format for disease prevalence analysis
long_data <- selected_data_75 %>%
  pivot_longer(cols = all_of(diseases),
               names_to = "Disease",
               values_to = "Status")

# Convert Status to binary (1 = Yes/1, 0 = No/0)
long_data <- long_data %>%
  mutate(Status = ifelse(Status == "Yes" | Status == "1" | Status == 1, 1, 0))

# Ensure Subgroup is a factor with correct levels
long_data$Subgroup <- factor(long_data$Subgroup, levels = c("Control", "C1", "C2", "C3", "C4", "C5"))

# Summarize disease prevalence: counts, proportions, percentages, and standard error
summary_data <- long_data %>%
  group_by(Subgroup, Disease) %>%
  summarise(
    Count = sum(Status),
    Total = n(),
    Proportion = Count / Total,
    SE = 1.96 * (sqrt((Proportion * (1 - Proportion)) / Total)),  # Standard error with 95% CI
    .groups = "drop"
  ) %>%
  mutate(Percentage = Proportion * 100)

# Define all subgroups for pairwise comparisons
all_subgroups <- c("Control", "C1", "C2", "C3", "C4", "C5")

# Calculate p-values for disease prevalence across subgroups
p_values <- long_data %>%
  filter(Subgroup %in% all_subgroups) %>%
  group_by(Disease) %>%
  group_modify(~ {
    # Remove missing values
    disease_data <- drop_na(.x)
    
    # Generate pairwise comparisons for all subgroup pairs
    pairwise_results <- map_dfr(combn(all_subgroups, 2, simplify = FALSE), function(subgroup_pair) {
      # Subset data for the current pair
      comparison_data <- disease_data %>%
        filter(Subgroup %in% subgroup_pair)
      
      # Create contingency table
      table_data <- table(comparison_data$Status, comparison_data$Subgroup)
      print(table_data)
      
      # Skip if table is incomplete (less than 2 columns or rows)
      if (ncol(table_data) < 2 || nrow(table_data) < 2) {
        return(tibble(Subgroup_1 = subgroup_pair[1], Subgroup_2 = subgroup_pair[2], p_value = NA))
      }
      
      # Perform Fisher’s test (for small counts) or Chi-square test
      p_value <- tryCatch({
        if (min(table_data) < 5) {
          fisher.test(table_data, simulate.p.value = TRUE, B = 1e5)$p.value
        } else {
          chisq.test(table_data)$p.value
        }
      }, error = function(e) NA)  # Return NA if test fails
      
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
  group_by(Disease) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))  # Apply Bonferroni correction

# Define custom labels for diseases
disease_labels <- c(
  "CKD1" = "Kidney Disease",
  "CVD1" = "Heart Disease",
  "Depression1" = "Depression",
  "Diabetes1" = "T2D",
  "Hypertension1" = "Hypertension",
  "Hyperlipidemia1" = "Hyperlipidemia",
  "Joint1" = "OA",
  "MetS1" = "MetS",
  "Migraine1" = "Migraine",
  "Obesity1" = "Obesity",
  "Sleep1" = "Sleep Apnea",
  "Fibrosis" = "Fibrosis",
  "Cirrhosis" = "Cirrhosis",
  "Transplant" = "Liver Transplant",
  "Hepatocellular" = "HCC"
)

# Order Disease factor by sorted custom labels
summary_data$Disease <- factor(summary_data$Disease, levels = names(sort(disease_labels)))

# Create grouped bar plot for disease prevalence
p <- ggplot(summary_data, aes(x = Disease, y = Percentage, fill = Subgroup)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, color = "gray40") +
  geom_errorbar(aes(ymin = Percentage - (SE * 100), ymax = Percentage + (SE * 100)),
                position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = c("#8B6914", "#B0C4DE", "#BBFFFF", 
                               "darkseagreen1", "lightsalmon", "plum1")) +
  scale_x_discrete(labels = disease_labels) +
  labs(x = "", y = "Percentage (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    legend.position = "none"
  )

# Save bar plot as PNG and PDF
ggsave("bio_disease_plot.png", plot = p, width = 10, height = 6, dpi = 400)
ggsave("bio_disease_plot.pdf", plot = p, width = 10, height = 6, dpi = 400)

# Export summary data and p-values to CSV
write.csv(summary_data, "disease_bio_propotion_results.csv", row.names = FALSE)
write.csv(p_values, "p_value_disease_bio_propotion_results.csv", row.names = FALSE)


# Conclusion
# ----------
# This script identifies latent clusters in MASLD biobank data using LCA, with results evaluated across
# multiple fit criteria. The code is designed for reproducibility and produces publication-ready outputs.


