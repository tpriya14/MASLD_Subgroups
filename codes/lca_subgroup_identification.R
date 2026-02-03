# ============================================================
# Latent Class Analysis Pipeline for MASLD Phenotyping
# ============================================================
# This script performs Latent Class Analysis (LCA) on MCB biobank data
# to identify clinically meaningful subgroups within the MASLD population.
#
# Overall Workflow:
#  1. Load and preprocess the dataset
#  2. Split MASLD cases into development and validation cohort
#  3. Fit LCA models on training data
#  4. Select and fit the final model
#  5. Assign latent class clusters
#  6. Generate summary tables and visualization datasets
#  7. Apply the LCA to the MCB and Tapestry whole dataset to use later as reference set
# ============================================================

# ------------------- Load Required Libraries -------------------
required_packages <- c(
  "data.table","dplyr","tidyr","purrr",
  "ggplot2","poLCA","MASS","table1","tidyLPA"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

run_blrt <- FALSE

# ------------------- Create Output Directories -------------------
dir.create("datasets", showWarnings = FALSE)
dir.create("datasets/lca_models", showWarnings = FALSE)
out_path <- "/datasets/lca_models/"

# ------------------- Load Main Dataset -------------------
final_all <- read.delim(
  "MCB_Data_Processed.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# ------------------- Prepare MASLD Cases for LCA -------------------
LCA_MASLD_data <- final_all[final_all$MASLD == 1, ]
selected_data_100 <- na.omit(LCA_MASLD_data)
LCA_data <- selected_data_100

character_columns <- sapply(LCA_data, is.character)
LCA_data[, character_columns] <- lapply(LCA_data[, character_columns], as.factor)

# ------------------- Split Data into Training and Test Sets -------------------
set.seed(5055)
test_index <- sample(1:nrow(selected_data_100), nrow(selected_data_100) * 0.25)
LCA_MASLD_data_75 <- selected_data_100[-test_index, ]
LCA_MASLD_data_25 <- selected_data_100[test_index, ]

# ------------------- Define Variables for LCA -------------------
vari_lca <- c(
  "PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", "MetS1",
  "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C",
  "Depression1", "Migraine1", "CKD1", "ALP_C"
)

# ------------------- Prepare LCA Input Datasets -------------------
dat_lca_75 <- LCA_MASLD_data_75[, vari_lca] %>%
  mutate(across(everything(), as.factor))

dat_lca_100 <- LCA_data
f_lca <- as.formula(paste0("cbind(", paste(vari_lca, collapse = ", "), ") ~ 1"))

# ------------------- Fit LCA Models on Training Data -------------------
n_start <- 1
n_end <- 10

for (i in n_start:n_end) {
  mi_start <- poLCA(
    f_lca, data = dat_lca_75, nclass = i,
    maxiter = 1000, nrep = 30,
    graph = FALSE, calc.se = TRUE, verbose = FALSE
  )
  
  vprobs.start <- poLCA.reorder(mi_start$probs.start,
                                order(mi_start$P, decreasing = TRUE))
  
  mi <- poLCA(
    f_lca, data = dat_lca_75, nclass = i,
    maxiter = 1000, probs.start = vprobs.start,
    graph = FALSE, calc.se = TRUE, verbose = FALSE
  )
  
  save(mi_start, mi, file = paste0(out_path, "lca_m", i, ".RData"))
}

save(dat_lca_75, file = paste0(out_path, "dat_lca_75.RData"))

# ------------------- Fit Final Selected Model (5 Classes) -------------------
optimal_k <- 5

lc_final <- poLCA(
  f_lca, data = dat_lca_75, nclass = optimal_k,
  maxiter = 1000, probs.start = vprobs.start,
  graph = FALSE, calc.se = TRUE, verbose = FALSE
)

save(lc_final, file = paste0(out_path, "lca_final_model.RData"))
LCA_MASLD_data_75$LatentClassCluster <- lc_final$predclass

# ------------------- Create Summary Table by Latent Class -------------------
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

exclude_cols <- c("PATIENT_ID", "MASLD")
table_columns <- LCA_MASLD_data_75[, !(names(LCA_MASLD_data_75) %in% exclude_cols)]

x <- table1(
  ~ . | LatentClassCluster,
  data = table_columns,
  overall = TRUE,
  render.categorical = my.render.cat,
  topclass = "Rtable1-zebra"
)

tab1_df <- as.data.frame(x)

write.table(tab1_df, "datasets/lca_summary_by_cluster.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------- Create Visualization Dataset (Development) -------------------
LCA_CTRL_data <- final_all[final_all$MASLD == 0, ]
LCA_CTRL_data_clean <- na.omit(LCA_CTRL_data)
LCA_CTRL_data_clean$LatentClassCluster <- 0

combined_data <- bind_rows(
  LCA_MASLD_data_75 %>% mutate(DatasetType = "Training"),
  LCA_CTRL_data_clean %>% mutate(DatasetType = "Control")
)

combined_data <- combined_data %>%
  mutate(Subgroup = case_when(
    LatentClassCluster == 0 ~ "Control",
    LatentClassCluster == 1 ~ "C1",
    LatentClassCluster == 2 ~ "C2",
    LatentClassCluster == 3 ~ "C3",
    LatentClassCluster == 4 ~ "C4",
    LatentClassCluster == 5 ~ "C5",
    TRUE ~ NA_character_
  ))

write.table(combined_data, "datasets/development_data_for_visualization.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------- Fit Final Model on Full (100%) Dataset -------------------
lc_final_100 <- poLCA(
  f_lca, data = dat_lca_100, nclass = optimal_k,
  maxiter = 1000, probs.start = vprobs.start,
  graph = FALSE, calc.se = TRUE, verbose = FALSE
)

save(lc_final_100, file = paste0(out_path, "lca_final_model_100.RData"))
LCA_MASLD_data_100 <- LCA_data
LCA_MASLD_data_100$LatentClassCluster <- lc_final_100$predclass

combined_data_test <- bind_rows(
  LCA_MASLD_data_100 %>% mutate(DatasetType = "Whole"),
  LCA_CTRL_data_clean %>% mutate(DatasetType = "Control")
)
combined_data_test <- combined_data_test %>%
  mutate(Subgroup = case_when(
    LatentClassCluster == 0 ~ "Control",
    LatentClassCluster == 1 ~ "C1",
    LatentClassCluster == 2 ~ "C2",
    LatentClassCluster == 3 ~ "C3",
    LatentClassCluster == 4 ~ "C4",
    LatentClassCluster == 5 ~ "C5",
    TRUE ~ NA_character_
  ))

# ------------------- Create Final Test Set with Clusters -------------------
LCA_MASLD_data_25 <- LCA_MASLD_data_100[test_index, ]
combined_data_test <- combined_data_test[test_index, ]
LCA_MASLD_data_25$Subgroup <- combined_data_test$LatentClassCluster

write.table(LCA_MASLD_data_100, "datasets/MASLD_100_with_clusters.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_data, "datasets/MASLD_75_with_clusters.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(LCA_MASLD_data_25, "datasets/MASLD_25_test_with_clusters.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_data_test, "datasets/test_data_for_visualization.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------- Apply Existing Model to TAPESTRY Dataset -------------------
tapestry_all <- read.delim("Tapestry_Data_Processed.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tapestry_clean <- na.omit(tapestry_all)
tapestry_clean <- tapestry_clean[tapestry_clean$MASLD == 1, ]
dat_lca_tapestry <- tapestry_clean[, vari_lca] %>% mutate(across(everything(), as.factor))

lc_tapestry <- poLCA(
  f_lca, data = dat_lca_tapestry, nclass = optimal_k,
  maxiter = 1000, probs.start = lc_final_100$probs.start,
  graph = FALSE, calc.se = TRUE, verbose = FALSE
)

tapestry_clean$LatentClassCluster <- lc_tapestry$predclass

LCA_CTRL_data_tap <- tapestry_all[tapestry_all$MASLD == 0, ]
LCA_CTRL_data_clean_tap <- na.omit(LCA_CTRL_data_tap)
LCA_CTRL_data_clean_tap$LatentClassCluster <- 0

combined_data_tap <- bind_rows(
  tapestry_clean %>% mutate(DatasetType = "Testing"),
  LCA_CTRL_data_clean_tap %>% mutate(DatasetType = "Control")
)

combined_data_tap <- combined_data_tap %>%
  mutate(Subgroup = case_when(
    LatentClassCluster == 0 ~ "Control",
    LatentClassCluster == 1 ~ "C1",
    LatentClassCluster == 2 ~ "C2",
    LatentClassCluster == 3 ~ "C3",
    LatentClassCluster == 4 ~ "C4",
    LatentClassCluster == 5 ~ "C5",
    TRUE ~ NA_character_
  ))

write.table(combined_data_tap, "datasets/Tapestry_with_clusters.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_data_tap, "datasets/tapestry_test_data_for_visualization.csv", sep = "\t", row.names = FALSE, quote = FALSE)
