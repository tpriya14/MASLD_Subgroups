# =============================================
# Core points-Based Membership Assignment
# =============================================
# This script assigns new patients to identified subgroups using DBSCAN-inspired approach.
# Core points are first identified. For each new patinet, the Euclidean distance 
# to the mean of each cluster's core points is calculated.
# Each sample is then assigned to the cluster with the smallest distance.
#
# Objective:
# - Assignment of MCB validation dataset
# - Assignment of independent validation Tapestry dataset
# - Model performance evaluation
# =============================================

# ------------------- Load Required Packages -------------------
required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", 
  "dbscan", "caret", "table1"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------- Load Training Data (75%) -------------------
selected_data_75 <- read.delim(
  "/datasets/MASLD_75_with_clusters.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# ------------------- Load and Prepare Test Data (25%) -------------------
test_data_25 <- read.delim(
  "/datasets/MASLD_25_test_with_clusters.csv",
  sep = "\t", header = TRUE
)

test_data_mcb <- test_data_25[, vari_lca]
char_cols <- sapply(test_data_mcb, is.character)
test_data_mcb[, char_cols] <- lapply(test_data_mcb[, char_cols], as.factor)

# ------------------- Load and Prepare Tapestry Cohort -------------------
tapestry_data <- read.delim(
  "/datasets/Tapestry_with_clusters.csv",
  sep = "\t", header = TRUE
)

tapestry_masld <- tapestry_data[tapestry_data$MASLD == 1, ]
final_all_tapestry <- tapestry_masld[, vari_lca]

char_cols_tap <- sapply(final_all_tapestry, is.character)
final_all_tapestry[, char_cols_tap] <- lapply(final_all_tapestry[, char_cols_tap], as.factor)

# ------------------- Define LCA Variables -------------------
vari_lca <- c(
  "LatentClassCluster", "PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", 
  "MetS1", "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C", 
  "Depression1", "Migraine1", "CKD1", "ALP_C"
)

# ------------------- Prepare Data with Dummy Variables -------------------
prepare_data_with_dummies <- function(data, vari_lca) {

  jacard_data <- data[, vari_lca]

  categorical_vars <- jacard_data %>% dplyr::select(-c(1))

  character_columns <- sapply(jacard_data, is.character)
  jacard_data[, character_columns] <- lapply(jacard_data[, character_columns], as.factor)

  dummy_categorical <- model.matrix(~ . - 1, data = categorical_vars)

  continuous_vars <- jacard_data %>% dplyr::select(c(1))
  combined_data <- cbind(continuous_vars, dummy_categorical)

  combined_data$LatentClassCluster <- jacard_data$LatentClassCluster

  return(combined_data)
}

train_data <- prepare_data_with_dummies(selected_data_75, vari_lca)

test_data_mcb_prep <- prepare_data_with_dummies(test_data_mcb, vari_lca)

test_data_tap_prep <- prepare_data_with_dummies(final_all_tapestry, vari_lca)

# ------------------- Align Columns Between Train and Test -------------------
align_test_with_train <- function(test_data, train_data) {

  train_cols <- colnames(train_data)
  test_cols <- colnames(test_data)

  missing_cols <- setdiff(train_cols, test_cols)
  if(length(missing_cols) > 0) {
    for(col in missing_cols) {
      test_data[[col]] <- 0
    }
  }

  extra_cols <- setdiff(test_cols, train_cols)
  if(length(extra_cols) > 0) {
    test_data <- test_data[, !(colnames(test_data) %in% extra_cols)]
  }

  test_data <- test_data[, train_cols]

  return(test_data)
}

test_data_mcb_prep <- align_test_with_train(test_data_mcb_prep, train_data)

test_data_tap_prep <- align_test_with_train(test_data_tap_prep, train_data)

# ------------------- DBSCAN Core Points Functions -------------------
get_core_points <- function(train_data, eps, minPts) {
  distances <- as.matrix(dist(train_data))
  neighbors <- apply(distances, 1, function(x) sum(x <= eps))
  core_points <- which(neighbors >= minPts)
  return(core_points)
}

kNNdistances <- function(data, k) {
  distances <- as.matrix(dist(data))
  kNN_distances <- apply(distances, 1, function(x) sort(x)[k + 1])
  return(sort(kNN_distances))
}

find_elbow_point <- function(distances) {
  n <- length(distances)
  second_derivative <- numeric(n)
  for (i in 2:(n-1)) {
    second_derivative[i] <- (distances[i+1] - 2 * distances[i] + distances[i-1])
  }
  elbow_point <- which.max(second_derivative)
  return(elbow_point)
}

# ------------------- Determine Optimal Eps for Each Cluster -------------------
k <- 13
minPts <- 13

num_clusters <- 5
cluster_eps <- numeric(num_clusters)

dir.create("/results/knn_plots", showWarnings = FALSE, recursive = TRUE)

for (i in 1:num_clusters) {

  train_cluster <- train_data[train_data$LatentClassCluster == i, ]

  if (nrow(train_cluster) < k + 1) {
    next
  }

  knn_distances <- kNNdistances(train_cluster[, -1], k)
  elbow_point <- find_elbow_point(knn_distances)
  eps_auto <- knn_distances[elbow_point]
  eps_final <- eps_auto
  cluster_eps[i] <- eps_final

  png(sprintf("validation_results/knn_plots/cluster_%d_knn_distances.png", i), width = 600, height = 400)

  plot(knn_distances, type = "l",
       main = sprintf("Cluster %d: k-NN Distance Plot (k=%d)", i, k),
       xlab = "Points (sorted by distance)",
       ylab = sprintf("%d-NN Distance", k))

  abline(h = eps_final, lty = 2)

  dev.off()
}

# ------------------- Compute Core Points for Each Cluster -------------------
cluster_core_points <- list()
cluster_avg <- list()

for (i in 1:num_clusters) {

  train_cluster <- train_data[train_data$LatentClassCluster == i, -1]

  if (nrow(train_cluster) == 0) next

  core_points <- get_core_points(train_cluster, cluster_eps[i], minPts)

  if (length(core_points) == 0) next

  cluster_core_points[[i]] <- core_points
  cluster_avg[[i]] <- colMeans(train_cluster[core_points, , drop = FALSE])
}

# ------------------- Prediction Function Using Core Points -------------------
predict_dbscan_core <- function(new_data, cluster_avg) {

  n_samples <- nrow(new_data)
  num_clusters <- length(cluster_avg)

  predicted_clusters <- integer(n_samples)

  for (j in 1:n_samples) {

    min_dist <- Inf
    cluster <- NA

    for (i in 1:num_clusters) {

      if (is.null(cluster_avg[[i]])) next

      avg <- cluster_avg[[i]]
      test_sample <- as.numeric(new_data[j, -c(1, ncol(new_data))])

      distance <- sqrt(sum((avg - test_sample)^2))

      if (distance < min_dist) {
        min_dist <- distance
        cluster <- i
      }
    }

    predicted_clusters[j] <- cluster
  }

  return(predicted_clusters)
}

# ------------------- Assign MCB Test Set -------------------
mcb_predictions <- predict_dbscan_core(test_data_mcb_prep, cluster_avg)

test_data_mcb$LatentClassCluster_DBSCAN <- mcb_predictions
test_data_mcb$Subgroup_DBSCAN <- paste0("C", mcb_predictions)

test_data_25$LatentClassCluster_DBSCAN <- mcb_predictions
test_data_25$Subgroup_DBSCAN <- paste0("C", mcb_predictions)

if("LatentClassCluster" %in% names(test_data_mcb)) {

  actual_mcb <- factor(test_data_mcb$LatentClassCluster, levels = 1:5)
  predicted_mcb <- factor(test_data_mcb$LatentClassCluster_DBSCAN, levels = 1:5)

  conf_matrix_mcb <- confusionMatrix(predicted_mcb, actual_mcb)
}

# ------------------- Assign Tapestry Cohort -------------------
tapestry_predictions <- predict_dbscan_core(test_data_tap_prep, cluster_avg)

final_all_tapestry$LatentClassCluster_DBSCAN <- tapestry_predictions
final_all_tapestry$Subgroup_DBSCAN <- paste0("C", tapestry_predictions)

tapestry_data$LatentClassCluster_DBSCAN <- tapestry_predictions
tapestry_data$Subgroup_DBSCAN <- paste0("C", tapestry_predictions)

if("LatentClassCluster" %in% names(test_data_tap_prep)) {

  actual_tap <- factor(test_data_tap_prep$LatentClassCluster, levels = 1:5)
  predicted_tap <- factor(tapestry_predictions, levels = 1:5)

  conf_matrix_tap <- confusionMatrix(predicted_tap, actual_tap)
}

# ------------------- Create Summary Tables -------------------
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

test_data_mcb_summary <- test_data_mcb[, !names(test_data_mcb) %in% c("PATIENT_ID", "MASLD", "Subgroup")]

mcb_table <- table1(~ . | LatentClassCluster_DBSCAN,
                    data = test_data_mcb_summary,
                    overall = TRUE,
                    render.categorical = my.render.cat,
                    topclass = "Rtable1-zebra")

mcb_table_df <- as.data.frame(mcb_table)

tapestry_summary <- final_all_tapestry[, !names(final_all_tapestry) %in% c("PATIENT_ID", "MASLD")]

tapestry_table <- table1(~ . | LatentClassCluster_DBSCAN,
                         data = tapestry_summary,
                         overall = TRUE,
                         render.categorical = my.render.cat,
                         topclass = "Rtable1-zebra")

tapestry_table_df <- as.data.frame(tapestry_table)

# ------------------- Save Results -------------------
write.table(
  test_data_25,
  "datasets/mcb_test_dbscan_assignments.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  mcb_table_df,
  "datasets/mcb_test_dbscan_summary.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
