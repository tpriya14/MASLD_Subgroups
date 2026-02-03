# =============================================
# Probability-Based Membership Assignment
# =============================================
# This script assign a new patient to latent classes based on their posterior probabilities, 
# which were calculated using the estimated class-conditional response probabilities 
# for the selected indicator variables derived from development cohort. 
#
# Objective:
# - Assignment of MCB validation dataset
# - Assignment of independent validation Tapestry dataset
# - Model performance evaluation
# =============================================

# ------------------- Load Required Packages -------------------
required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr",
  "poLCA", "table1", "caret", "ggplot2"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

setwd("/")

# ------------------- Custom LCA Functions -------------------

poLCA.vectorize <- function(probs) {
  classes <- nrow(probs[[1]])
  vecprobs <- unlist(lapply(probs, t))
  numChoices <- sapply(probs, ncol)
  return(list(vecprobs = vecprobs, numChoices = numChoices, classes = classes))
}

poLCA.updatePrior <- function(coeff, x, nclass) {
  prior <- matrix(0, nrow = nrow(x), ncol = nclass)
  for(k in 1:nclass) {
    prior[, k] <- ifelse(k == 1, 1, exp(x %*% coeff[[k]]))
  }
  prior <- prior / rowSums(prior)
  return(prior)
}

poLCA.posterior <- function(lc, y, x = NULL) {

  if (is.vector(y) || any(dim(y) == 1)) y <- matrix(y, nrow = 1)
  y[is.na(y)] <- 0

  if (!is.null(x)) {
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    x <- cbind(1, x)
  }

  prior <- if (!is.null(x) && ncol(lc$x) > 1) {
    poLCA.updatePrior(lc$coeff, x, length(lc$P))
  } else {
    matrix(lc$P, nrow = nrow(y), ncol = length(lc$P), byrow = TRUE)
  }

  ret <- poLCA.postClass.C(prior, poLCA.vectorize(lc$probs), y)
  return(ret)
}

poLCA.postClass.C <- function(prior, vp, y) {

  ret <- .C("postclass",
            as.double(t(prior)),
            as.double(vp$vecprobs),
            as.integer(t(y)),
            as.integer(length(vp$numChoices)),
            as.integer(dim(y)[1]),
            as.integer(vp$numChoices),
            as.integer(vp$classes),
            posterior = double(dim(y)[1] * vp$classes))

  ret$posterior <- matrix(ret$posterior, ncol = vp$classes, byrow = TRUE)

  return(ret$posterior)
}

predict_lc <- function(lc, new_data, formula) {

  mframe <- model.frame(formula, new_data, na.action = na.pass)
  y <- model.response(mframe)

  posterior <- poLCA.posterior(lc = lc, y = y)
  predicted_class <- apply(posterior, 1, which.max)

  return(list(
    class = predicted_class,
    posterior = posterior,
    max_prob = apply(posterior, 1, max)
  ))
}

# ------------------- Load Training Data and LCA Model -------------------

load("/datasets/lca_models/lca_final_model.RData")

training_data <- read.delim(
  "/datasets/MASLD_75_with_clusters.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

vari_lca <- c(
  "PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", "MetS1", 
  "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C", 
  "Depression1", "Migraine1", "CKD1", "ALP_C"
)

f_lca <- as.formula(paste0("cbind(", paste(vari_lca, collapse = ", "), ") ~ 1"))

# ------------------- Load and Prepare Test Data (25%) -------------------

test_data_25 <- read.delim(
  "/datasets/MASLD_25_test_with_clusters.csv",
  sep = "\t", header = TRUE
)

test_lca_data <- test_data_25[, vari_lca]

char_cols <- sapply(test_lca_data, is.character)
test_lca_data[, char_cols] <- lapply(test_lca_data[, char_cols], as.factor)

# ------------------- Load and Prepare Tapestry Cohort -------------------

tapestry_data <- read.delim(
  "/datasets/Tapestry_with_clusters.csv",
  sep = "\t", header = TRUE
)

tapestry_masld <- tapestry_data[tapestry_data$MASLD == 1, ]

tapestry_lca_data <- tapestry_masld[, vari_lca]

char_cols_tap <- sapply(tapestry_lca_data, is.character)
tapestry_lca_data[, char_cols_tap] <- lapply(tapestry_lca_data[, char_cols_tap], as.factor)

# ------------------- Predict Latent Classes for Test Set -------------------

test_predictions <- predict_lc(lc_final, test_lca_data, f_lca)

test_data_25 <- test_data_25 %>%
  mutate(
    LatentClassCluster_Predicted = test_predictions$class,
    PredictionProbability = test_predictions$max_prob,
    Subgroup_Predicted = case_when(
      LatentClassCluster_Predicted == 1 ~ "C1",
      LatentClassCluster_Predicted == 2 ~ "C2",
      LatentClassCluster_Predicted == 3 ~ "C3",
      LatentClassCluster_Predicted == 4 ~ "C4",
      LatentClassCluster_Predicted == 5 ~ "C5"
    )
  )

# ------------------- Predict Latent Classes for Tapestry Cohort -------------------

tapestry_predictions <- predict_lc(lc_final, tapestry_lca_data, f_lca)

tapestry_masld <- tapestry_masld %>%
  mutate(
    LatentClassCluster_Predicted = tapestry_predictions$class,
    PredictionProbability = tapestry_predictions$max_prob,
    Subgroup_Predicted = case_when(
      LatentClassCluster_Predicted == 1 ~ "C1",
      LatentClassCluster_Predicted == 2 ~ "C2",
      LatentClassCluster_Predicted == 3 ~ "C3",
      LatentClassCluster_Predicted == 4 ~ "C4",
      LatentClassCluster_Predicted == 5 ~ "C5"
    )
  )

# ------------------- Create Summary Tables -------------------

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

test_summary_cols <- test_data_25[, !names(test_data_25) %in% 
                                    c("PATIENT_ID", "MASLD", "Subgroup", 
                                      "LatentClassCluster", "Subgroup_Predicted")]

test_table <- table1(
  ~ . | LatentClassCluster_Predicted,
  data = test_summary_cols,
  overall = TRUE,
  render.categorical = my.render.cat,
  topclass = "Rtable1-zebra"
)

test_table_df <- as.data.frame(test_table)

tapestry_summary_cols <- tapestry_masld[, !names(tapestry_masld) %in% 
                                          c("PATIENT_ID", "MASLD", "LatentClassCluster")]

tapestry_table <- table1(
  ~ . | Subgroup,
  data = tapestry_summary_cols,
  overall = TRUE,
  render.categorical = my.render.cat,
  topclass = "Rtable1-zebra"
)

tapestry_table_df <- as.data.frame(tapestry_table)

# ------------------- Model Performance Evaluation -------------------

dir.create("validation_results", showWarnings = FALSE)

if("LatentClassCluster" %in% names(test_data_25)) {

  actual_test <- factor(test_data_25$LatentClassCluster, levels = 1:5, labels = paste0("C", 1:5))
  predicted_test <- factor(test_data_25$Subgroup_Predicted, levels = paste0("C", 1:5))

  conf_matrix_test <- confusionMatrix(predicted_test, actual_test)
}

if("LatentClassCluster" %in% names(tapestry_masld)) {

  actual_test <- factor(tapestry_masld$LatentClassCluster, levels = 1:5, labels = paste0("C", 1:5))
  predicted_test <- factor(tapestry_masld$Subgroup_Predicted, levels = paste0("C", 1:5))

  conf_matrix_test <- confusionMatrix(predicted_test, actual_test)
}

# ------------------- Save Results -------------------

write.table(
  test_data_25,
  "/datasets/prob_validation_test_set_with_predictions.csv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  test_table_df,
  "/datasets/prob_validation_set_summary_by_cluster.csv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  tapestry_masld,
  "/datasets/prob_validation_tapestry_with_predictions.csv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  tapestry_table_df,
  "/datasets/prob_validation_tapestry_summary_by_cluster.csv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
