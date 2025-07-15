# Title: Probability-Based Membership Assignment
# Description: This script performs Probability-Based Membership Assignment using clinical data to assign patients to subgroups.

all_pkgs <- unique(c(
  "data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer", # Overall
  "poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA", # LCA-specific
  "olsrr", "dplyr", "tidyr", "purrr", "tidyfit", "table1", # Additional
  "RcppArmadillo", "BranchGLM", "MASS", "caret", "mlr3misc" # Additional
))

# Load packages
sapply(all_pkgs, require, character.only = TRUE, quietly = TRUE)


# ------------------- Set Working Directory -------------------
setwd("Journal_Final_Pic")

# ------------------- Custom LCA Functions -------------------
# Function to vectorize LCA probabilities
poLCA.vectorize <- function(probs) {
  classes <- nrow(probs[[1]])
  vecprobs <- unlist(lapply(probs, t))
  numChoices <- sapply(probs, ncol)
  return(list(vecprobs = vecprobs, numChoices = numChoices, classes = classes))
}

# Function to compute posterior probabilities for LCA
poLCA.posterior <- function(lc, y, x = NULL) {
  if (is.vector(y) | any(dim(y) == 1)) {
    y <- matrix(y, nrow = 1)
  }
  y[is.na(y)] <- 0
  
  if (!is.null(x)) {
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    x <- cbind(1, x)
  }
  
  if ((ncol(lc$x) > 1) & (!is.null(x))) {
    prior <- poLCA.updatePrior(lc$coeff, x, length(lc$P))
  } else {
    prior <- matrix(lc$P, nrow = nrow(y), ncol = length(lc$P), byrow = TRUE)
  }
  ret <- poLCA.postClass.C(prior, poLCA.vectorize(lc$probs), y)
  return(ret)
}

# Function to call C routine for posterior class probabilities
poLCA.postClass.C <- function(prior, vp, y) {
  ret <- .C("postclass",
            as.double(t(prior)),
            as.double(vp$vecprobs),
            as.integer(t(y)),
            as.integer(length(vp$numChoices)),
            as.integer(dim(y)[1]),
            as.integer(vp$numChoices),
            as.integer(vp$classes),
            posterior = double(dim(y)[1] * vp$classes)
  )
  ret$posterior <- matrix(ret$posterior, ncol = vp$classes, byrow = TRUE)
  return(ret$posterior)
}

# Function to predict latent class memberships
predict_lc <- function(lc, new_data) {
  mframe <- model.frame(f_lca, new_data)
  y <- model.response(mframe)
  n_classes <- 5
  posterior <- poLCA.posterior(lc = lc, y = y)
  pr_cl <- apply(posterior, 1, which.max)
  return(pr_cl)
}

# ------------------- Data Loading and Preprocessing -------------------
# Load datasets
selected_data_75 <- read.delim("T75_5_Bio_After_LCA_T1533_clinical_info_06_10.csv", 
                               sep = "\t", header = TRUE)
selected_data_25 <- read.delim("T25_5_Bio_After_LCA_T510_clinical_info_06_11.csv", 
                               sep = "\t", header = TRUE)

# Inspect training data
nrow(selected_data_75)
colnames(selected_data_75)

# ------------------- LCA Model Setup -------------------
# Define variables for LCA
vari_lca <- c("PATIENT_GENDER_NAME", "Obesity1", "Hyperlipidemia1", "MetS1", 
              "Hypertension1", "ALT_C", "AST_C", "BMI_C", "HDL_C", 
              "Depression1", "Migraine1", "CKD1", "ALP_C")

# Create LCA formula
f_lca <- eval(parse(text = paste0("as.formula(cbind(", 
                                  paste(vari_lca, collapse = ", "), ") ~ 1)")))

# Prepare test data for LCA
test_data <- selected_data_25[, vari_lca]
character_columns <- sapply(test_data, is.character)
test_data[, character_columns] <- lapply(test_data[, character_columns], as.factor)
test_data$LatentClassCluster <- selected_data_25$LatentClassCluster
test_data$LatentClassCluster1 <- NA
test_data_1 <- test_data[, -c(ncol(test_data) - 1, ncol(test_data))]

# Load pre-trained LCA model
load("lc5.RData")

# ------------------- LCA Prediction -------------------
# Predict latent class memberships
predictions <- predict_lc(lc5, test_data_1)
test_data$LatentClassCluster1 <- predictions
selected_data_25$LatentClassCluster1 <- test_data$LatentClassCluster1

# Assign subgroup labels
selected_data_25 <- selected_data_25 %>% 
  mutate(Subgroup = case_when(
    LatentClassCluster1 == 1 ~ "C1",
    LatentClassCluster1 == 2 ~ "C2",
    LatentClassCluster1 == 3 ~ "C3",
    LatentClassCluster1 == 4 ~ "C4",
    LatentClassCluster1 == 5 ~ "C5"
  ))

# ------------------- Summary Tables and Evaluation -------------------
# Custom rendering function for categorical variables
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%0.1f", PCT))))
}

# Create summary table by LatentClassCluster1
x <- table1(~ . | LatentClassCluster1,
            data = selected_data_25[, -c(1:4)], overall = TRUE, 
            render.categorical = my.render.cat, topclass = "Rtable1-zebra")
tab1_df <- as.data.frame(x)

# Evaluate predictions with confusion matrix
actual <- factor(selected_data_25$LatentClassCluster)
predicted <- factor(selected_data_25$LatentClassCluster1)
conf_matrix <- confusionMatrix(actual, predicted)

# ------------------- Save Outputs -------------------
# Save updated test data
write.table(selected_data_25, 
            file = "T25_5_Bio_After_LCA_T510_clinical_info_10_14.csv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save confusion matrix
save(conf_matrix, file = "confusion_matrix_Bio_25_LCA_10_14.RData")

# Save summary table
write.table(tab1_df, file = "NTap_After_LCA_matched_Bio.csv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ------------------- End of Script -------------------
