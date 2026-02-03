# =============================================
# Centroid-Based Membership Assignment
# =============================================
# This script assigns patients to the identified subgroups from the development cohort 
# using latent class analysis.
#
# Tasks:
# - Assignment of MCB validation dataset
# - Assignment of independent validation Tapestry dataset
# - Model performance evaluation
#
# =============================================

# Load Required Packages
# ----------------------
pkgs_overall <- c("data.table", "reshape2", "ggplot2", "dplyr", "tidyr", "caret")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)

cat("========================================\n")
cat("Centroid-Based Cluster Assignment\n")
cat("========================================\n\n")

# ------------------- Load Training Data (75%) -------------------
selected_data_75 <- read.delim(
  "/datasets/MASLD_75_with_clusters.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
nrow(selected_data_75)
# ------------------------------------------------------------
# 4. Load and Prepare Test Data (25%)
# ------------------------------------------------------------
cat("Step 2: Loading 25% test set\n")
cat("----------------------------\n")

test_data_25 <- read.delim(
  "/datasets/MASLD_25_test_with_clusters.csv",
  sep = "\t", header = TRUE
)
cat("✓ Test set loaded:", nrow(test_data_25), "samples\n")

# Prepare for LCA prediction
test_data_mcb <- test_data_25[, vari_lca]
char_cols <- sapply(test_data_mcb, is.character)
test_data_mcb[, char_cols] <- lapply(test_data_mcb[, char_cols], as.factor)
cat("✓ Test data prepared for LCA prediction\n\n")

# ------------------------------------------------------------
# 5. Load and Prepare Tapestry Cohort
# ------------------------------------------------------------
cat("Step 3: Loading Tapestry validation cohort\n")
cat("------------------------------------------\n")

tapestry_data <- read.delim(
  "/datasets/Tapestry_with_clusters.csv",
  sep = "\t", header = TRUE
)

tapestry_masld <- tapestry_data[tapestry_data$MASLD == 1, ]
final_all_tapestry <- tapestry_masld[, vari_lca]
char_cols_tap <- sapply(final_all_tapestry, is.character)
final_all_tapestry[, char_cols_tap] <- lapply(final_all_tapestry[, char_cols_tap], as.factor)
cat("✓ Tapestry cohort prepared:", nrow(tapestry_masld), "MASLD samples\n\n")

# ------------------- Define LCA Variables -------------------
vari_lca <- c(
  "LatentClassCluster", "PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", 
  "MetS1", "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C", 
  "Depression1", "Migraine1", "CKD1", "ALP_C"
)

# ------------------- Prepare Training Data -------------------

dat_lca_75 <- selected_data_75[, vari_lca] %>%
  filter(MASLD==1)
train_data <- dat_lca_75
jacard_data_bio_75 <- train_data
categorical_vars <- jacard_data_bio_75 %>% dplyr::select(-c(1))
character_columns <- sapply(jacard_data_bio_75, is.character)
jacard_data_bio_75[, character_columns] <- lapply(
  jacard_data_bio_75[, character_columns], 
  as.factor
)
dummy_categorical <- model.matrix(~ . - 1, data = categorical_vars)
continuous_vars <- jacard_data_bio_75 %>% dplyr::select(c(1))
combined_data_bio_75 <- cbind(continuous_vars, dummy_categorical)
combined_data_bio_75$LatentClassCluster <- jacard_data_bio_75$LatentClassCluster
# ------------------- Function to Prepare Test Data -------------------
prepare_test_data <- function(test_dataset, vari_lca) {
  test_data_subset <- test_dataset[, vari_lca]
  jacard_data_test <- test_data_subset
  categorical_vars <- jacard_data_test %>% dplyr::select(-c(1))
  character_columns <- sapply(jacard_data_test, is.character)
  jacard_data_test[, character_columns] <- lapply(
    jacard_data_test[, character_columns], 
    as.factor
  )
  dummy_categorical <- model.matrix(~ . - 1, data = categorical_vars)
  continuous_vars <- jacard_data_test %>% dplyr::select(c(1))
  combined_data_test <- cbind(continuous_vars, dummy_categorical)
  combined_data_test$LatentClassCluster <- jacard_data_test$LatentClassCluster
  
  return(combined_data_test)
}

# ------------------- Prepare MCB Test Data -------------------
combined_data_mcb <- prepare_test_data(test_data_mcb, vari_lca)
colnames(test_data_mcb)
train_cols <- colnames(combined_data_bio_75)
test_cols <- colnames(combined_data_mcb)
missing_cols <- setdiff(train_cols, test_cols)
extra_cols <- setdiff(test_cols, train_cols)
combined_data_mcb <- combined_data_mcb[, train_cols]

# ------------------- Prepare Tapestry Test Data -------------------

combined_data_tapestry <- prepare_test_data(final_all_tapestry, vari_lca)
test_cols_tap <- colnames(combined_data_tapestry)
missing_cols_tap <- setdiff(train_cols, test_cols_tap)
extra_cols_tap <- setdiff(test_cols_tap, train_cols)
combined_data_tapestry <- combined_data_tapestry[, train_cols]

# ------------------- Distance Calculation Function -------------------
calculate_KMeans_similarity <- function(test_row_idx, cluster_num, train_data, test_data) {
  centroids <- colMeans(train_data[train_data$LatentClassCluster == cluster_num, -c(1)])
  test_sample <- test_data[test_row_idx, -c(1,22)]
  distance <- sqrt(sum((centroids - test_sample)^2))
  return(distance)
}

# ------------------- Assign MCB Test Set -------------------
cat("========================================\n")
cat("Assigning MCB Test Set (25%)\n")
cat("========================================\n\n")

train_data <- combined_data_bio_75
test_data <- combined_data_mcb
test_data$LatentClassCluster1 <- NA
print(colnames(train_data))
print(colnames(test_data))
num_clusters <- 5
differences <- numeric(length = nrow(test_data))
pb <- txtProgressBar(min = 0, max = nrow(test_data), style = 3)

for (i in 1:nrow(test_data)) {
  distances <- numeric(length = num_clusters)
  min_dist <- Inf
  assigned_cluster <- NA
  for (j in 1:num_clusters) {
    jaccard_val <- calculate_KMeans_similarity(i, j, train_data, subset(test_data, select = -LatentClassCluster1))
    distances[j] <- jaccard_val
    if (jaccard_val < min_dist) {
      min_dist <- jaccard_val
      assigned_cluster <- j
    }
  }
  sorted_arr <- sort(distances)
  differences[i] <- sorted_arr[2] - sorted_arr[1]
  test_data[i, ]$LatentClassCluster1 <- assigned_cluster
  setTxtProgressBar(pb, i)
}
close(pb)

test_data_mcb$LatentClassCluster_Centroid <- test_data$LatentClassCluster1
test_data_mcb$Subgroup_Centroid <- paste0("C", test_data$LatentClassCluster1)
print(table(test_data_mcb$Subgroup_Centroid))
cat("\n")

# ------------------- Evaluate MCB Performance -------------------
actual <- factor(test_data$LatentClassCluster, levels = 1:5)
predicted <- factor(test_data$LatentClassCluster1, levels = 1:5)
conf_matrix_mcb <- confusionMatrix(actual, predicted)
cat("\nConfusion Matrix - MCB Test Set:\n")
print(conf_matrix_mcb)
# ------------------- Assign Tapestry Cohort -------------------

cat("========================================\n")
cat("Assigning Tapestry Cohort\n")
cat("========================================\n\n")

test_data_tap <- combined_data_tapestry
test_data_tap$LatentClassCluster1 <- NA

differences_tap <- numeric(length = nrow(test_data_tap))
pb <- txtProgressBar(min = 0, max = nrow(test_data_tap), style = 3)

for (i in 1:nrow(test_data_tap)) {
  distances <- numeric(length = num_clusters)
  min_dist <- Inf
  assigned_cluster <- NA
  for (j in 1:num_clusters) {
    jaccard_val <- calculate_KMeans_similarity(i, j, train_data, subset(test_data_tap, select = -LatentClassCluster1))
    distances[j] <- jaccard_val
    if (jaccard_val < min_dist) {
      min_dist <- jaccard_val
      assigned_cluster <- j
    }
  }
  sorted_arr <- sort(distances)
  differences_tap[i] <- sorted_arr[2] - sorted_arr[1]
  test_data_tap[i, ]$LatentClassCluster1 <- assigned_cluster
  setTxtProgressBar(pb, i)
}
close(pb)
final_all_tapestry$LatentClassCluster_Centroid <- test_data_tap$LatentClassCluster1
final_all_tapestry$Subgroup_Centroid <- paste0("C", test_data_tap$LatentClassCluster1)

actual_tap <- factor(test_data_tap$LatentClassCluster, levels = 1:5)
predicted_tap <- factor(test_data_tap$LatentClassCluster1, levels = 1:5)
conf_matrix_tap <- confusionMatrix(actual_tap, predicted_tap)
cat("Confusion Matrix - Tapestry Cohort:\n")
print(conf_matrix_tap)



# ------------------- Save Results -------------------
cat("========================================\n")
cat("Saving Results\n")
cat("========================================\n\n")

test_data_25$LatentClassCluster_Centroid <- test_data$LatentClassCluster1
test_data_25$Subgroup_Centroid <- paste0("C", test_data$LatentClassCluster1)

# Save MCB test set results
write.table(
  test_data_25,
  "datasets/centroid_mcb_test_centroid_assignments.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

nrow(selected_data_75)
colnames(selected_data_75)
colnames(test_data_25)
test_data_25$LatentClassCluster <- test_data_25$Subgroup_Predicted
combined_mcb <- rbind(selected_data_75[,-c(50)], test_data_25[, -c(50:52)])

# Save MCB test set results
write.table(
  combined_mcb,
  "datasets/centroid_mcb_1A_1B_test_centroid_assignments.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# Save Tapestry results
tapestry_data$LatentClassCluster_Centroid <- final_all_tapestry$LatentClassCluster_Centroid
tapestry_data$Subgroup_Centroid <- paste0("C", final_all_tapestry$Subgroup_Centroid)

write.table(
  tapestry_data,
  "datasets/centroid_tapestry_centroid_assignments.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
