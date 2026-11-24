
# ============================================================
# LIGHTGBM MULTICLASS CLASSIFICATION
# 8 GROUPS | 24 SAMPLES | 3 REPLICATES/GROUP
# 10 REPEATS OF STRATIFIED 3-FOLD CROSS-VALIDATION
# ============================================================


# 1. INSTALL AND LOAD PACKAGES
# ============================================================

if (!requireNamespace("lightgbm", quietly = TRUE)) {
  install.packages("lightgbm")
}

if (!requireNamespace("lightgbm", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  remotes::install_github(
    "lightgbm-org/LightGBM",
    subdir = "R-package"
  )
}

packages <- c("lightgbm", "ggplot2", "dplyr",
              "tidyr", "pheatmap", "pROC")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(lightgbm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

set.seed(123)


# 2. IMPORT YOUR CSV DATA
# ============================================================

data <- read.csv(
  file.choose(),
  header = FALSE,
  sep = ",",
  fill = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Imported data dimensions:\n")
print(dim(data))

# First row: sample names
sample_names <- as.character(unlist(data[1, -1]))

# Second row: class labels (0-7)
group_labels <- as.character(unlist(data[2, -1]))

# Remaining rows: protein IDs and abundances
protein_ids <- as.character(unlist(data[-c(1, 2), 1]))

abundance_data <- data[-c(1, 2), -1, drop = FALSE]

abundance_data[] <- lapply(
  abundance_data,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

abundance_matrix <- as.matrix(abundance_data)

# Remove rows without protein IDs
valid_ids <- !is.na(protein_ids) & trimws(protein_ids) != ""

protein_ids <- protein_ids[valid_ids]
abundance_matrix <- abundance_matrix[valid_ids, , drop = FALSE]

# Remove completely missing and all-zero proteins
keep <- rowSums(!is.na(abundance_matrix)) > 0 &
  rowSums(abundance_matrix, na.rm = TRUE) != 0

protein_ids <- protein_ids[keep]
abundance_matrix <- abundance_matrix[keep, , drop = FALSE]

# Ensure unique protein IDs
protein_ids <- make.unique(protein_ids)

# Convert to samples x proteins
X <- t(abundance_matrix)

rownames(X) <- sample_names
colnames(X) <- protein_ids

cat("\nSamples x proteins:\n")
print(dim(X))


# 3. DEFINE THE EIGHT CLASSES
# ============================================================

class_levels <- as.character(0:7)

group <- factor(
  group_labels,
  levels = class_levels
)

if (any(is.na(group))) {
  stop("Invalid group labels. Check the second row of your CSV.")
}

if (nrow(X) != length(group)) {
  stop("Sample number and group-label number do not match.")
}

cat("\nGroup distribution:\n")
print(table(group))

if (any(table(group) != 3)) {
  stop("Expected three replicates in each of the eight groups.")
}


# 4. OUTPUT FOLDER
# ============================================================

output_dir <- "LightGBM_Results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 5. PREPROCESSING FUNCTION
# ============================================================

# Tree models do not require scaling.
# Genuine zero values are preserved.
# Missing values are imputed using training-set medians.

preprocess_data <- function(train_x, test_x) {
  
  train_x <- as.matrix(train_x)
  test_x <- as.matrix(test_x)
  
  # Remove proteins completely missing in training data
  valid <- colSums(!is.na(train_x)) > 0
  
  train_x <- train_x[, valid, drop = FALSE]
  test_x <- test_x[, valid, drop = FALSE]
  
  # Impute missing values using training medians
  for (j in seq_len(ncol(train_x))) {
    
    med <- median(train_x[, j], na.rm = TRUE)
    
    train_x[is.na(train_x[, j]), j] <- med
    test_x[is.na(test_x[, j]), j] <- med
  }
  
  # Remove zero-variance proteins
  feature_sd <- apply(train_x, 2, sd)
  
  keep <- is.finite(feature_sd) & feature_sd > 0
  
  train_x <- train_x[, keep, drop = FALSE]
  test_x <- test_x[, keep, drop = FALSE]
  
  if (ncol(train_x) == 0) {
    stop("No variable proteins remain after preprocessing.")
  }
  
  # LightGBM requires numeric matrices
  storage.mode(train_x) <- "double"
  storage.mode(test_x) <- "double"
  
  list(
    train = train_x,
    test = test_x
  )
}


# 6. STRATIFIED FOLD GENERATION
# ============================================================

make_stratified_folds <- function(y, k = 3) {
  
  folds <- integer(length(y))
  
  for (lev in levels(y)) {
    
    ids <- which(y == lev)
    ids <- sample(ids)
    
    folds[ids] <- rep(seq_len(k), length.out = length(ids))
  }
  
  folds
}


# 7. PERFORMANCE METRICS
# ============================================================

calculate_metrics <- function(actual, predicted, levels) {
  
  actual <- factor(actual, levels = levels)
  predicted <- factor(predicted, levels = levels)
  
  cm <- table(
    Actual = actual,
    Predicted = predicted
  )
  
  accuracy <- sum(diag(cm)) / sum(cm)
  
  precision <- recall <- f1 <- numeric(length(levels))
  
  for (i in seq_along(levels)) {
    
    tp <- cm[i, i]
    fp <- sum(cm[, i]) - tp
    fn <- sum(cm[i, ]) - tp
    
    precision[i] <- ifelse(
      tp + fp == 0, 0, tp / (tp + fp)
    )
    
    recall[i] <- ifelse(
      tp + fn == 0, 0, tp / (tp + fn)
    )
    
    f1[i] <- ifelse(
      precision[i] + recall[i] == 0,
      0,
      2 * precision[i] * recall[i] /
        (precision[i] + recall[i])
    )
  }
  
  list(
    Accuracy = accuracy,
    Precision = mean(precision),
    Recall = mean(recall),
    F1 = mean(f1),
    Confusion = cm
  )
}


# 7B. MULTICLASS ROC-AUC (ONE-VS-REST)
# AUC is calculated separately for each class using out-of-fold
# probabilities pooled across the 3 folds of each repeat.
# Macro-AUC is the unweighted mean of the 8 class-specific AUCs.

calculate_ovr_auc <- function(actual, probability_matrix, class_levels) {
  actual <- as.character(actual)
  probability_matrix <- as.matrix(probability_matrix)

  if (ncol(probability_matrix) != length(class_levels)) {
    stop("Probability matrix must contain one column per class.")
  }
  colnames(probability_matrix) <- class_levels

  class_auc <- setNames(rep(NA_real_, length(class_levels)),
                        paste0("AUC_Class_", class_levels))

  for (lev in class_levels) {
    binary_actual <- factor(
      ifelse(actual == lev, "Target", "Other"),
      levels = c("Other", "Target")
    )

    roc_obj <- tryCatch(
      pROC::roc(
        response = binary_actual,
        predictor = probability_matrix[, lev],
        levels = c("Other", "Target"),
        direction = "<",
        quiet = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(roc_obj)) {
      class_auc[paste0("AUC_Class_", lev)] <-
        as.numeric(pROC::auc(roc_obj))
    }
  }

  macro_auc <- if (all(is.na(class_auc))) NA_real_ else
    mean(class_auc, na.rm = TRUE)

  list(Macro_AUC = macro_auc, Class_AUC = class_auc)
}


# 8. LIGHTGBM PARAMETERS
# ============================================================

# Conservative parameters for a very small dataset.
# Avoid using a large number of leaves.

params <- list(
  objective = "multiclass",
  metric = "multi_logloss",
  num_class = 8,
  
  learning_rate = 0.05,
  num_leaves = 3,
  max_depth = 2,
  min_data_in_leaf = 2,
  
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 1,
  
  lambda_l1 = 1,
  lambda_l2 = 5,
  
  verbosity = -1,
  num_threads = 2,
  seed = 123
)

nrounds <- 60


# 9. INITIALIZE RESULTS
# ============================================================

prediction_results <- list()
repeat_metrics <- list()
importance_results <- list()

prediction_counter <- 1
importance_counter <- 1

n_repeats <- 10
n_folds <- 3


# 10. REPEATED STRATIFIED CROSS-VALIDATION
# Macro one-vs-rest AUC is calculated from held-out class probabilities.
# ============================================================

for (repeat_id in seq_len(n_repeats)) {
  
  cat("\n================================\n")
  cat("REPEAT:", repeat_id, "\n")
  cat("================================\n")
  
  fold_assignment <- make_stratified_folds(
    group,
    k = n_folds
  )
  
  repeat_predictions <- list()
  
  for (fold in seq_len(n_folds)) {
    
    cat("Repeat", repeat_id, "- Fold", fold, "\n")
    
    train_index <- which(fold_assignment != fold)
    test_index <- which(fold_assignment == fold)
    
    x_train_raw <- X[train_index, , drop = FALSE]
    x_test_raw <- X[test_index, , drop = FALSE]
    
    y_train <- group[train_index]
    y_test <- group[test_index]
    
    # Preprocessing
    processed <- preprocess_data(
      x_train_raw,
      x_test_raw
    )
    
    train_matrix <- processed$train
    test_matrix <- processed$test
    
    # LightGBM labels must be integers from 0 to 7
    train_label <- as.integer(as.character(y_train))
    test_label <- as.integer(as.character(y_test))
    
    # Training dataset
    dtrain <- lgb.Dataset(
      data = train_matrix,
      label = train_label
    )
    
    # Train model
    model <- lgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds
    )
    
    # Predict class probabilities
    prob <- predict(
      model,
      test_matrix
    )
    
    prob <- as.matrix(prob)
    
    if (ncol(prob) != 8) {
      stop("LightGBM did not return eight class probabilities.")
    }
    
    # Predicted class
    predicted <- max.col(
      prob,
      ties.method = "first"
    ) - 1
    
    predicted <- as.character(predicted)
    
    # Save predictions
    probability_df <- as.data.frame(prob)
    colnames(probability_df) <- paste0("Prob_Class_", class_levels)

    fold_predictions <- data.frame(
      Repeat = repeat_id,
      Fold = fold,
      Sample = rownames(test_matrix),
      Actual = as.character(test_label),
      Predicted = predicted,
      probability_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    prediction_results[[prediction_counter]] <-
      fold_predictions
    
    repeat_predictions[[fold]] <- fold_predictions
    
    prediction_counter <- prediction_counter + 1
    
    # Extract feature importance
    importance <- lgb.importance(
      model = model,
      percentage = FALSE
    )
    
    if (!is.null(importance) && nrow(importance) > 0) {
      
      importance$Repeat <- repeat_id
      importance$Fold <- fold
      
      importance_results[[importance_counter]] <-
        importance
      
      importance_counter <- importance_counter + 1
    }
  }
  
  # Repeat-wise performance
  repeat_df <- bind_rows(repeat_predictions)
  
  metrics <- calculate_metrics(
    repeat_df$Actual,
    repeat_df$Predicted,
    class_levels
  )
  
  repeat_probability_matrix <- as.matrix(
    repeat_df[, paste0("Prob_Class_", class_levels), drop = FALSE]
  )

  auc_metrics <- calculate_ovr_auc(
    actual = repeat_df$Actual,
    probability_matrix = repeat_probability_matrix,
    class_levels = class_levels
  )

  repeat_metrics[[repeat_id]] <- data.frame(
    Repeat = repeat_id,
    Accuracy = metrics$Accuracy,
    Precision = metrics$Precision,
    Recall = metrics$Recall,
    F1 = metrics$F1,
    Macro_AUC = auc_metrics$Macro_AUC,
    as.list(auc_metrics$Class_AUC),
    check.names = FALSE
  )
}


# 11. COMBINE PREDICTIONS
# ============================================================

all_predictions <- bind_rows(prediction_results)

cat("\nPrediction dimensions:\n")
print(dim(all_predictions))

cat("\nPrediction columns:\n")
print(names(all_predictions))

if (nrow(all_predictions) == 0) {
  stop("No predictions generated.")
}

write.csv(
  all_predictions,
  file.path(output_dir, "LightGBM_All_Predictions.csv"),
  row.names = FALSE
)


# 12. PERFORMANCE SUMMARY
# ============================================================

metrics_df <- bind_rows(repeat_metrics)

write.csv(
  metrics_df,
  file.path(output_dir, "LightGBM_Repeat_Metrics.csv"),
  row.names = FALSE
)

performance_summary <- metrics_df %>%
  summarise(
    Accuracy_Mean = mean(Accuracy),
    Accuracy_SD = sd(Accuracy),
    
    Precision_Mean = mean(Precision),
    Precision_SD = sd(Precision),
    
    Recall_Mean = mean(Recall),
    Recall_SD = sd(Recall),
    
    F1_Mean = mean(F1),
    F1_SD = sd(F1),

    Macro_AUC_Mean = mean(Macro_AUC, na.rm = TRUE),
    Macro_AUC_SD = sd(Macro_AUC, na.rm = TRUE)
  )

cat("\nLIGHTGBM PERFORMANCE:\n")
print(performance_summary)

write.csv(
  performance_summary,
  file.path(output_dir, "LightGBM_Performance_Summary.csv"),
  row.names = FALSE
)


# 13. POOLED CONFUSION MATRIX
# ============================================================

pooled_cm <- table(
  Actual = factor(
    all_predictions$Actual,
    levels = class_levels
  ),
  Predicted = factor(
    all_predictions$Predicted,
    levels = class_levels
  )
)

write.csv(
  as.data.frame.matrix(pooled_cm),
  file.path(output_dir, "LightGBM_Confusion_Matrix_Counts.csv")
)

cm_percent <- prop.table(
  pooled_cm,
  margin = 1
) * 100

write.csv(
  as.data.frame.matrix(round(cm_percent, 2)),
  file.path(output_dir, "LightGBM_Confusion_Matrix_Percent.csv")
)

cat("\nPOOLED CONFUSION MATRIX:\n")
print(pooled_cm)

png(
  file.path(output_dir, "LightGBM_Confusion_Matrix.png"),
  width = 2200,
  height = 1800,
  res = 250
)

pheatmap(
  cm_percent,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  main = "LightGBM Classification\nRow-normalized Confusion Matrix (%)",
  angle_col = 45
)

dev.off()


# 14. PERFORMANCE PLOT
# ============================================================

performance_long <- metrics_df %>%
  pivot_longer(
    cols = c(Accuracy, Precision, Recall, F1, Macro_AUC),
    names_to = "Metric",
    values_to = "Value"
  )

performance_plot <- ggplot(
  performance_long,
  aes(x = Metric, y = Value)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    title = "LightGBM Classification Performance",
    x = "Performance metric",
    y = "Score"
  )

ggsave(
  file.path(output_dir, "LightGBM_Performance.png"),
  performance_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# 15. COMBINE FEATURE IMPORTANCE FROM CV FOLDS
# ============================================================

if (length(importance_results) > 0) {
  
  importance_df <- bind_rows(importance_results)
  
  write.csv(
    importance_df,
    file.path(output_dir, "LightGBM_Foldwise_Importance.csv"),
    row.names = FALSE
  )
  
  # Aggregate importance across all 30 folds
  # Gain = contribution to model improvement
  # Frequency = number of times a protein is used in splits
  
  importance_summary <- importance_df %>%
    group_by(Feature) %>%
    summarise(
      Mean_Gain = mean(Gain, na.rm = TRUE),
      Total_Gain = sum(Gain, na.rm = TRUE),
      Mean_Split = mean(Frequency, na.rm = TRUE),
      Total_Split = sum(Frequency, na.rm = TRUE),
      Folds_Selected = n_distinct(
        paste(Repeat, Fold)
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Gain))
  
  importance_summary$Selection_Frequency_Percent <-
    importance_summary$Folds_Selected /
    (n_repeats * n_folds) * 100
  
  write.csv(
    importance_summary,
    file.path(output_dir, "LightGBM_Protein_Importance_CV.csv"),
    row.names = FALSE
  )
  
  cat("\nTOP 20 LIGHTGBM PROTEINS FROM CV:\n")
  print(head(importance_summary, 20))
  
  top20 <- head(importance_summary, 20)
  
  importance_plot <- ggplot(
    top20,
    aes(
      x = reorder(Feature, Mean_Gain),
      y = Mean_Gain
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      title = "Top 20 Proteins - LightGBM",
      x = "Protein ID",
      y = "Mean Gain across CV folds"
    )
  
  ggsave(
    file.path(output_dir, "LightGBM_Top20_Proteins_CV.png"),
    importance_plot,
    width = 8,
    height = 7,
    dpi = 300
  )
  
} else {
  
  cat("\nNo feature importance was returned from CV models.\n")
}


# 16. FINAL MODEL USING ALL 24 SAMPLES
# ============================================================

# This model is for exploratory feature ranking only.
# Its importance is not independent validation.

final_processed <- preprocess_data(X, X)

final_X <- final_processed$train

final_label <- as.integer(as.character(group))

final_train <- lgb.Dataset(
  data = final_X,
  label = final_label
)

final_model <- lgb.train(
  params = params,
  data = final_train,
  nrounds = nrounds
)

saveRDS(
  final_model,
  file.path(output_dir, "Final_LightGBM_Model.rds")
)

final_importance <- lgb.importance(
  model = final_model,
  percentage = TRUE
)

write.csv(
  final_importance,
  file.path(output_dir, "LightGBM_Final_Model_Importance.csv"),
  row.names = FALSE
)

cat("\nTOP 20 FINAL MODEL PROTEINS:\n")
print(head(final_importance, 20))


# 17. FINISH
# ============================================================

cat("\n====================================\n")
cat("LIGHTGBM ANALYSIS COMPLETED\n")
cat("====================================\n")

cat("\nTotal cross-validation predictions:\n")
print(nrow(all_predictions))

cat("\nPerformance summary:\n")
print(performance_summary)

cat("\nOutput folder:\n")
print(normalizePath(output_dir))

cat("\nAnalysis completed successfully.\n")

