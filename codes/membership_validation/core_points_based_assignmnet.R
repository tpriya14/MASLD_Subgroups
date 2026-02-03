# =============================================
# Core points-Based Membership Assignment
# =============================================
# This script assigns new patients to identified subgroups using DBSCAN-inspired approach.
# Core points are first identified. For each new patinet, the Euclidean distance 
# to the mean of each cluster's core points is calculated.
# Each sample is then assigned to the cluster with the smallest distance.
#
# Tasks:
# - Assignment of MCB validation dataset
# - Assignment of independent validation Tapestry dataset
# - Model performance evaluation
# =============================================

# =============================================
# Load Required Packages
# =============================================
required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", 
  "dbscan", "caret", "table1"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all required libraries
invisible(lapply(required_packages, library, character.only = TRUE))

cat("========================================\n")
cat("DBSCAN Core Points Assignment\n")
cat("========================================\n\n")

# =============================================
# Step 1: Load Training Data (75%)
# =============================================
selected_data_75 <- read.delim(
  "/datasets/MASLD_75_with_clusters.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# =============================================
# Step 2: Load and Prepare Test Data (25%)
# =============================================
cat("Step 2: Loading 25% test set\n")
cat("----------------------------\n")

test_data_25 <- read.delim(
  "/datasets/MASLD_25_test_with_clusters.csv",
  sep = "\t", header = TRUE
)
cat("✓ Test set loaded:", nrow(test_data_25), "samples\n")

# Prepare test data for LCA prediction
test_data_mcb <- test_data_25[, vari_lca]
char_cols <- sapply(test_data_mcb, is.character)
test_data_mcb[, char_cols] <- lapply(test_data_mcb[, char_cols], as.factor)
cat("✓ Test data prepared for LCA prediction\n\n")

# =============================================
# Step 3: Load and Prepare Tapestry Cohort
# =============================================
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

# =============================================
# Step 4: Define LCA Variables
# =============================================
vari_lca <- c(
  "LatentClassCluster", "PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", 
  "MetS1", "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C", 
  "Depression1", "Migraine1", "CKD1", "ALP_C"
)

# =============================================
# Step 5: Prepare Data with Dummy Variables
# =============================================
cat("Step 4: Creating dummy variables\n")
cat("---------------------------------\n")

prepare_data_with_dummies <- function(data, vari_lca) {
  # Select LCA variables
  jacard_data <- data[, vari_lca]
  
  # Separate cluster labels from features
  categorical_vars <- jacard_data %>% dplyr::select(-c(1))
  
  # Convert character columns to factors
  character_columns <- sapply(jacard_data, is.character)
  jacard_data[, character_columns] <- lapply(jacard_data[, character_columns], as.factor)
  
  # Create dummy variables for categorical columns
  dummy_categorical <- model.matrix(~ . - 1, data = categorical_vars)
  
  # Combine with continuous variables
  continuous_vars <- jacard_data %>% dplyr::select(c(1))
  combined_data <- cbind(continuous_vars, dummy_categorical)
  combined_data$LatentClassCluster <- jacard_data$LatentClassCluster
  
  return(combined_data)
}

# Prepare training and test datasets
train_data <- prepare_data_with_dummies(selected_data_75, vari_lca)
cat("✓ Training data prepared:", ncol(train_data) - 1, "features\n")

test_data_mcb_prep <- prepare_data_with_dummies(test_data_mcb, vari_lca)
cat("✓ MCB test data prepared\n")

test_data_tap_prep <- prepare_data_with_dummies(final_all_tapestry, vari_lca)
cat("✓ Tapestry data prepared\n\n")

# =============================================
# Step 6: Align Columns Between Train and Test
# =============================================
cat("Step 5: Aligning columns between train and test\n")
cat("-----------------------------------------------\n")

align_test_with_train <- function(test_data, train_data) {
  train_cols <- colnames(train_data)
  test_cols <- colnames(test_data)
  
  # Add missing columns as 0
  missing_cols <- setdiff(train_cols, test_cols)
  if(length(missing_cols) > 0) {
    for(col in missing_cols) {
      test_data[[col]] <- 0
    }
  }
  
  # Remove extra columns
  extra_cols <- setdiff(test_cols, train_cols)
  if(length(extra_cols) > 0) {
    test_data <- test_data[, !(colnames(test_data) %in% extra_cols)]
  }
  
  # Reorder columns to match training data
  test_data <- test_data[, train_cols]
  
  return(test_data)
}

test_data_mcb_prep <- align_test_with_train(test_data_mcb_prep, train_data)
cat("✓ MCB test data aligned\n")

test_data_tap_prep <- align_test_with_train(test_data_tap_prep, train_data)
cat("✓ Tapestry data aligned\n\n")

# =============================================
# Step 7: Define DBSCAN Core Points Functions
# =============================================
cat("Step 6: Setting up DBSCAN core points detection\n")
cat("-----------------------------------------------\n")

# Identify core points based on DBSCAN concept
get_core_points <- function(train_data, eps, minPts) {
  distances <- as.matrix(dist(train_data))
  neighbors <- apply(distances, 1, function(x) sum(x <= eps))
  core_points <- which(neighbors >= minPts)
  return(core_points)
}

# Compute k-nearest neighbors distances
kNNdistances <- function(data, k) {
  distances <- as.matrix(dist(data))
  kNN_distances <- apply(distances, 1, function(x) sort(x)[k + 1])
  return(sort(kNN_distances))
}

# Find elbow point for optimal eps
find_elbow_point <- function(distances) {
  n <- length(distances)
  second_derivative <- numeric(n)
  for (i in 2:(n-1)) {
    second_derivative[i] <- (distances[i+1] - 2 * distances[i] + distances[i-1])
  }
  elbow_point <- which.max(second_derivative)
  return(elbow_point)
}

cat("✓ DBSCAN functions defined\n\n")

# =============================================
# Step 8: Determine Optimal Eps for Each Cluster
# =============================================
cat("Step 7: Determining optimal eps for each cluster\n")
cat("-------------------------------------------------\n")

k <- 13  # Number of neighbors for kNN
minPts <- 13  # Minimum points for DBSCAN core points

num_clusters <- 5
cluster_eps <- numeric(num_clusters)

dir.create("/results/knn_plots", 
           showWarnings = FALSE, recursive = TRUE)

for (i in 1:num_clusters) {
  train_cluster <- train_data[train_data$LatentClassCluster == i, ]
  
  if (nrow(train_cluster) < k + 1) {
    cat("  Cluster", i, ": Not enough samples\n")
    next
  }
  
  knn_distances <- kNNdistances(train_cluster[, -1], k)
  elbow_point <- find_elbow_point(knn_distances)
  eps_auto <- knn_distances[elbow_point]
  
  # Manual adjustment for specific clusters
  if (i == 1 || i == 3) {
    eps_final <- 1.5
  } else if (i == 2) {
    eps_final <- 1.8
  } else if (i == 4 || i == 5) {
    eps_final <- 2.0
  } else {
    eps_final <- eps_auto
  }
  
  cluster_eps[i] <- eps_final
  cat(sprintf("  Cluster %d: eps = %.2f (auto: %.2f)\n", i, eps_final, eps_auto))
  
  # Save kNN distance plot
  png(sprintf("validation_results/knn_plots/cluster_%d_knn_distances.png", i), width = 600, height = 400)
  plot(knn_distances, type = "l", 
       main = sprintf("Cluster %d: k-NN Distance Plot (k=%d)", i, k),
       xlab = "Points (sorted by distance)", 
       ylab = sprintf("%d-NN Distance", k),
       col = "blue", lwd = 2)
  abline(h = eps_final, col = "red", lty = 2, lwd = 2)
  legend("topleft", legend = c("k-NN distances", sprintf("eps = %.2f", eps_final)),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2)
  dev.off()
}

cat("\n✓ Optimal eps values determined\n\n")

# =============================================
# Step 9: Compute Core Points for Each Cluster
# =============================================
cat("Step 8: Computing core points for each cluster\n")
cat("-----------------------------------------------\n")

cluster_core_points <- list()
cluster_avg <- list()

for (i in 1:num_clusters) {
  train_cluster <- train_data[train_data$LatentClassCluster == i, -1]
  
  if (nrow(train_cluster) == 0) next
  
  core_points <- get_core_points(train_cluster, cluster_eps[i], minPts)
  
  if (length(core_points) == 0) {
    cat(sprintf("  Cluster %d: No core points found\n", i))
    next
  }
  
  cluster_core_points[[i]] <- core_points
  cluster_avg[[i]] <- colMeans(train_cluster[core_points, , drop = FALSE])
  
  cat(sprintf("  Cluster %d: %d core points (%.1f%% of cluster)\n", 
              i, length(core_points), 100 * length(core_points) / nrow(train_cluster)))
}

cat("\n✓ Core points computed for all clusters\n\n")

# =============================================
# Step 10: Prediction Function Using Core Points
# =============================================
predict_dbscan_core <- function(new_data, cluster_avg) {
  n_samples <- nrow(new_data)
  num_clusters <- length(cluster_avg)
  
  predicted_clusters <- integer(n_samples)
  distances_matrix <- matrix(NA, nrow = n_samples, ncol = num_clusters)
  
  cat("  Assigning samples...\n")
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  
  for (j in 1:n_samples) {
    distances_c <- numeric(num_clusters)
    min_dist <- Inf
    cluster <- NA
    
    for (i in 1:num_clusters) {
      if (is.null(cluster_avg[[i]])) {
        distances_c[i] <- Inf
        next
      }
      
      avg <- cluster_avg[[i]]
      test_sample <- as.numeric(new_data[j, -c(1, ncol(new_data))])
      distance <- sqrt(sum((avg - test_sample)^2))
      distances_c[i] <- distance
      
      if (distance < min_dist) {
        min_dist <- distance
        cluster <- i
      }
    }
    
    predicted_clusters[j] <- cluster
    distances_matrix[j, ] <- distances_c
    
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  return(list(
    clusters = predicted_clusters,
    distances = distances_matrix
  ))
}

# =============================================
# Step 11: Assign MCB Test Set
# =============================================
cat("Step 11: Assigning MCB test set (25%)\n")
cat("=====================================\n")

mcb_predictions <- predict_dbscan_core(test_data_mcb_prep, cluster_avg)

# Add predictions to original dataset
test_data_mcb$LatentClassCluster_DBSCAN <- mcb_predictions$clusters
test_data_mcb$Subgroup_DBSCAN <- paste0("C", mcb_predictions$clusters)
test_data_25$LatentClassCluster_DBSCAN <- mcb_predictions$clusters
test_data_25$Subgroup_DBSCAN <- paste0("C", mcb_predictions$clusters)

cat("\n✓ MCB test set assigned\n")
cat("\nClass distribution:\n")
print(table(test_data_mcb$Subgroup_DBSCAN))

# Evaluate performance if labels exist
if("LatentClassCluster" %in% names(test_data_mcb)) {
  actual_mcb <- factor(test_data_mcb$LatentClassCluster, levels = 1:5)
  predicted_mcb <- factor(test_data_mcb$LatentClassCluster_DBSCAN, levels = 1:5)
  
  conf_matrix_mcb <- confusionMatrix(predicted_mcb, actual_mcb)
  
  cat("\nConfusion Matrix - MCB Test Set:\n")
  print(conf_matrix_mcb$table)
  cat("\nPerformance:\n")
  cat("  Accuracy:", round(conf_matrix_mcb$overall["Accuracy"], 3), "\n")
  cat("  Kappa:", round(conf_matrix_mcb$overall["Kappa"], 3), "\n\n")
}

# =============================================
# Step 12: Assign Tapestry Cohort
# =============================================
cat("Step 12: Assigning Tapestry cohort\n")
cat("==================================\n")

tapestry_predictions <- predict_dbscan_core(test_data_tap_prep, cluster_avg)

final_all_tapestry$LatentClassCluster_DBSCAN <- tapestry_predictions$clusters
final_all_tapestry$Subgroup_DBSCAN <- paste0("C", tapestry_predictions$clusters)
tapestry_data$LatentClassCluster_DBSCAN <- tapestry_predictions$clusters
tapestry_data$Subgroup_DBSCAN <- paste0("C", tapestry_predictions$clusters)

cat("\n✓ Tapestry cohort assigned\n")
cat("\nClass distribution:\n")
print(table(final_all_tapestry$Subgroup_DBSCAN))

# Evaluate if true labels exist
if("LatentClassCluster" %in% names(test_data_tap_prep) && 
   !all(is.na(test_data_tap_prep$LatentClassCluster))) {
  actual_tap <- factor(test_data_tap_prep$LatentClassCluster, levels = 1:5)
  predicted_tap <- factor(tapestry_predictions$clusters, levels = 1:5)
  
  conf_matrix_tap <- confusionMatrix(predicted_tap, actual_tap)
  
  cat("\nConfusion Matrix - Tapestry:\n")
  print(conf_matrix_tap$table)
  cat("\nPerformance:\n")
  cat("  Accuracy:", round(conf_matrix_tap$overall["Accuracy"], 3), "\n")
  cat("  Kappa:", round(conf_matrix_tap$overall["Kappa"], 3), "\n\n")
}

# =============================================
# Step 13: Create Summary Tables
# =============================================
cat("Step 13: Creating summary tables\n")
cat("---------------------------------\n")

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

# MCB Test Set Summary
test_data_mcb_summary <- test_data_mcb[, !names(test_data_mcb) %in% c("PATIENT_ID", "MASLD", "Subgroup")]
mcb_table <- table1(~ . | LatentClassCluster_DBSCAN, data = test_data_mcb_summary,
                    overall = TRUE, render.categorical = my.render.cat, topclass = "Rtable1-zebra")
mcb_table_df <- as.data.frame(mcb_table)
cat("✓ MCB summary table created\n")

# Tapestry Summary
tapestry_summary <- final_all_tapestry[, !names(final_all_tapestry) %in% c("PATIENT_ID", "MASLD")]
tapestry_table <- table1(~ . | LatentClassCluster_DBSCAN, data = tapestry_summary,
                         overall = TRUE, render.categorical = my.render.cat, topclass = "Rtable1-zebra")
tapestry_table_df <- as.data.frame(tapestry_table)
cat("✓ Tapestry summary table created\n\n")

# =============================================
# Step 14: Save Results
# =============================================
cat("Step 14: Saving results\n")
cat("-----------------------\n")

# Save MCB assignments and summary
write.table(test_data_25, "datasets/mcb_test_dbscan_assignments.csv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ MCB assignments saved\n")
write.table(mcb_table_df, "datasets/mcb_test_dbscan_summary.csv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ MCB summary table saved\n")

# Save Tapestry assignments and summary
write.table(tapestry_data, "datasets/tapestry_dbscan_assignments.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tapestry_table_df, "datasets/tapestry_dbscan_summary.csv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ Tapestry assignments saved\n")

                                               
