
############################################################
# MULTICLASS LINEAR SVM ANALYSIS
# Dataset: 8 experimental groups, 3 replicates/group
# Method: Repeated stratified 3-fold cross-validation
# Repeats: 10
# Outputs: overall and class-wise metrics, balanced accuracy, confusion matrix
############################################################

# 1. INSTALL AND LOAD PACKAGES
packages <- c("e1071", "ggplot2", "dplyr", "tidyr", "pROC")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)

set.seed(123)

# 2. SELECT AND IMPORT YOUR CSV FILE

file_path <- file.choose()

data <- read.csv(
  file_path,
  header = FALSE,
  sep = ",",
  fill = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Imported rows:", nrow(data), "\n")
cat("Imported columns:", ncol(data), "\n")

if (ncol(data) != 25) {
  stop("Expected 25 columns: 1 protein ID + 24 samples.")
}

# 3. EXTRACT SAMPLE NAMES AND GROUP LABELS

sample_names <- as.character(unlist(data[1, -1]))
sample_names <- trimws(sample_names)

labels <- as.character(unlist(data[2, -1]))
labels <- trimws(labels)

group_names <- c(
  "Fe3O4_Positive_In_vivo",
  "Fe3O4_Positive_In_vitro",
  "Fe3O4_Negative_In_vivo",
  "Fe3O4_Negative_In_vitro",
  "Green_Hard",
  "Chem_Hard",
  "Green_Soft",
  "Chem_Soft"
)

label_map <- setNames(group_names, as.character(0:7))

if (length(labels) != 24) {
  stop("The dataset must contain 24 sample labels.")
}

if (!all(labels %in% names(label_map))) {
  stop("Unexpected group labels found in row 2.")
}

group <- factor(
  unname(label_map[labels]),
  levels = group_names
)

if (any(table(group) != 3)) {
  stop("Expected exactly 3 replicates in each of the 8 groups.")
}

cat("\nSample distribution:\n")
print(table(group))

# 4. EXTRACT PROTEIN ABUNDANCE MATRIX

protein_data <- data[-c(1, 2), , drop = FALSE]

protein_ids <- trimws(as.character(protein_data[[1]]))

# Convert abundance columns to numeric
abundance <- as.matrix(
  data.frame(
    lapply(protein_data[-1], function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    }),
    check.names = FALSE
  )
)

storage.mode(abundance) <- "numeric"

rownames(abundance) <- make.unique(protein_ids)

# Remove rows without protein IDs
valid_ids <- !is.na(protein_ids) & protein_ids != ""
abundance <- abundance[valid_ids, , drop = FALSE]

# Remove proteins that are zero in every sample
# Zeros in individual samples are retained
nonzero_proteins <- rowSums(abundance != 0, na.rm = TRUE) > 0
abundance <- abundance[nonzero_proteins, , drop = FALSE]

# Remove proteins that contain only missing values
valid_proteins <- rowSums(!is.na(abundance)) > 0
abundance <- abundance[valid_proteins, , drop = FALSE]

cat("\nProteins retained:", nrow(abundance), "\n")

# Samples must be rows and proteins must be columns
X <- t(abundance)

rownames(X) <- make.unique(sample_names)
colnames(X) <- make.unique(rownames(abundance))

storage.mode(X) <- "numeric"

cat("Samples:", nrow(X), "\n")
cat("Protein features:", ncol(X), "\n")

# 5. FUNCTION FOR TRAINING-BASED PREPROCESSING
# Imputation and scaling are learned from training samples only.
# Zero abundance values are not converted to missing values.

prepare_data <- function(train_x, test_x) {
  
  train_x <- as.matrix(train_x)
  test_x <- as.matrix(test_x)

  # Convert NA, NaN, Inf and -Inf to missing values.
  # Imputation parameters are still estimated from training data only.
  train_x[!is.finite(train_x)] <- NA_real_
  test_x[!is.finite(test_x)] <- NA_real_
  
  # Remove features that are entirely missing in training data
  keep <- colSums(!is.na(train_x)) > 0
  
  train_x <- train_x[, keep, drop = FALSE]
  test_x <- test_x[, keep, drop = FALSE]
  
  if (ncol(train_x) == 0) {
    stop("No usable protein features in a training fold.")
  }
  
  # Impute missing values using training medians
  medians <- apply(train_x, 2, median, na.rm = TRUE)
  
  for (j in seq_len(ncol(train_x))) {
    
    train_x[is.na(train_x[, j]), j] <- medians[j]
    test_x[is.na(test_x[, j]), j] <- medians[j]
    
  }
  
  # Remove zero-variance features in training data
  train_sd <- apply(train_x, 2, sd)
  
  keep_sd <- is.finite(train_sd) & train_sd > 0
  
  train_x <- train_x[, keep_sd, drop = FALSE]
  test_x <- test_x[, keep_sd, drop = FALSE]
  
  if (ncol(train_x) == 0) {
    stop("No variable protein features remain in a training fold.")
  }
  
  # Scale using training means and standard deviations
  train_means <- colMeans(train_x)
  train_sd <- apply(train_x, 2, sd)
  
  train_scaled <- scale(
    train_x,
    center = train_means,
    scale = train_sd
  )
  
  test_scaled <- scale(
    test_x,
    center = train_means,
    scale = train_sd
  )
  
  return(list(
    train = as.matrix(train_scaled),
    test = as.matrix(test_scaled)
  ))
}

# 6. METRIC CALCULATION FUNCTION

calculate_metrics <- function(actual, predicted, group_names) {
  
  cm <- table(
    factor(actual, levels = group_names),
    factor(predicted, levels = group_names)
  )
  
  precision <- numeric(length(group_names))
  recall <- numeric(length(group_names))
  f1 <- numeric(length(group_names))
  
  for (i in seq_along(group_names)) {
    
    TP <- cm[i, i]
    FP <- sum(cm[, i]) - TP
    FN <- sum(cm[i, ]) - TP
    
    precision[i] <- if ((TP + FP) > 0) {
      TP / (TP + FP)
    } else {
      0
    }
    
    recall[i] <- if ((TP + FN) > 0) {
      TP / (TP + FN)
    } else {
      0
    }
    
    f1[i] <- if ((precision[i] + recall[i]) > 0) {
      2 * precision[i] * recall[i] /
        (precision[i] + recall[i])
    } else {
      0
    }
  }
  
  accuracy <- sum(diag(cm)) / sum(cm)
  balanced_accuracy <- mean(recall)
  
  return(c(
    Accuracy = accuracy,
    Balanced_Accuracy = balanced_accuracy,
    Macro_Precision = mean(precision),
    Macro_Recall = mean(recall),
    Macro_F1 = mean(f1)
  ))
}

# Calculate precision, recall and F1 separately for every class.
calculate_class_metrics <- function(actual, predicted, group_names) {

  cm <- table(
    factor(actual, levels = group_names),
    factor(predicted, levels = group_names)
  )

  class_results <- data.frame(
    Class = group_names,
    Precision = numeric(length(group_names)),
    Recall = numeric(length(group_names)),
    F1 = numeric(length(group_names)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(group_names)) {
    TP <- cm[i, i]
    FP <- sum(cm[, i]) - TP
    FN <- sum(cm[i, ]) - TP

    class_results$Precision[i] <- if ((TP + FP) > 0) TP / (TP + FP) else 0
    class_results$Recall[i] <- if ((TP + FN) > 0) TP / (TP + FN) else 0

    p <- class_results$Precision[i]
    r <- class_results$Recall[i]
    class_results$F1[i] <- if ((p + r) > 0) 2 * p * r / (p + r) else 0
  }

  class_results
}

# One-vs-rest ROC-AUC for each class and unweighted macro-AUC.
# Probabilities must be from held-out samples only.
calculate_ovr_auc <- function(actual, probability_matrix, group_names) {

  auc_values <- setNames(
    rep(NA_real_, length(group_names)),
    paste0("AUC_Class_", 0:(length(group_names) - 1))
  )

  if (is.null(probability_matrix) || nrow(probability_matrix) != length(actual)) {
    warning("Probability matrix is missing or does not match the number of held-out samples.")
    return(c(Macro_AUC = NA_real_, auc_values))
  }

  for (i in seq_along(group_names)) {
    target_class <- group_names[i]
    binary_actual <- as.integer(as.character(actual) == target_class)

    if (length(unique(binary_actual)) == 2 &&
        target_class %in% colnames(probability_matrix)) {

      roc_obj <- tryCatch(
        pROC::roc(
          response = binary_actual,
          predictor = as.numeric(probability_matrix[, target_class]),
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        ),
        error = function(e) NULL
      )

      if (!is.null(roc_obj)) {
        auc_values[i] <- as.numeric(pROC::auc(roc_obj))
      }
    }
  }

  valid_auc <- auc_values[is.finite(auc_values)]
  macro_auc <- if (length(valid_auc) > 0) mean(valid_auc) else NA_real_

  c(Macro_AUC = macro_auc, auc_values)
}

# 7. REPEATED STRATIFIED 3-FOLD CROSS-VALIDATION

n_folds <- 3
n_repeats <- 10

all_predictions <- data.frame()
repeat_metrics <- data.frame()
classwise_repeat_metrics <- data.frame()

for (repeat_id in seq_len(n_repeats)) {
  
  cat("\nRunning repeat:", repeat_id, "of", n_repeats, "\n")
  
  fold_id <- integer(length(group))
  
  # Assign one replicate from each group to each fold
  for (g in group_names) {
    
    group_indices <- which(group == g)
    
    shuffled_indices <- sample(group_indices)
    
    fold_id[shuffled_indices] <- rep(
      seq_len(n_folds),
      length.out = length(shuffled_indices)
    )
  }
  
  repeat_actual <- character()
  repeat_predicted <- character()
  repeat_probabilities <- NULL
  
  for (fold in seq_len(n_folds)) {
    
    test_indices <- which(fold_id == fold)
    train_indices <- which(fold_id != fold)
    
    train_x <- X[train_indices, , drop = FALSE]
    test_x <- X[test_indices, , drop = FALSE]
    
    train_y <- group[train_indices]
    test_y <- group[test_indices]
    
    # Preprocess using training data only
    processed <- prepare_data(train_x, test_x)
    
    # Linear SVM
    svm_model <- svm(
      x = processed$train,
      y = train_y,
      type = "C-classification",
      kernel = "linear",
      cost = 1,
      scale = FALSE,
      probability = TRUE
    )

    # Predict class labels and estimate class probabilities for held-out samples.
    predictions_raw <- predict(
      svm_model,
      newdata = processed$test,
      probability = TRUE
    )

    fold_probabilities <- attr(predictions_raw, "probabilities")

    if (is.null(fold_probabilities)) {
      stop("SVM probabilities were not returned. Check probability = TRUE.")
    }

    # Reorder probability columns to match group_names exactly.
    probability_matrix <- matrix(
      NA_real_,
      nrow = nrow(processed$test),
      ncol = length(group_names),
      dimnames = list(NULL, group_names)
    )

    common_classes <- intersect(colnames(fold_probabilities), group_names)
    probability_matrix[, common_classes] <-
      fold_probabilities[, common_classes, drop = FALSE]

    if (any(!is.finite(probability_matrix))) {
      stop("Missing or non-finite SVM class probabilities detected.")
    }

    # Append probabilities in the same order as repeat_actual.
    repeat_probabilities <- rbind(
      repeat_probabilities,
      probability_matrix
    )

    predictions <- factor(
      as.character(predictions_raw),
      levels = group_names
    )
    
    # Save predictions
    fold_results <- data.frame(
      Repeat = repeat_id,
      Fold = fold,
      Sample = rownames(X)[test_indices],
      Actual = as.character(test_y),
      Predicted = as.character(predictions)
    )
    
    all_predictions <- rbind(
      all_predictions,
      fold_results
    )
    
    repeat_actual <- c(
      repeat_actual,
      as.character(test_y)
    )
    
    repeat_predicted <- c(
      repeat_predicted,
      as.character(predictions)
    )
  }
  
  # Confirm that each held-out sample has one row of class probabilities.
  if (is.null(repeat_probabilities) ||
      nrow(repeat_probabilities) != length(repeat_actual)) {
    stop("Probability rows do not match held-out observations in this repeat.")
  }

  # Calculate metrics for this complete repeat
  metrics <- calculate_metrics(
    repeat_actual,
    repeat_predicted,
    group_names
  )
  
  auc_metrics <- calculate_ovr_auc(
    actual = repeat_actual,
    probability_matrix = repeat_probabilities,
    group_names = group_names
  )

  # Class-wise precision, recall, F1 and one-vs-rest AUC for this repeat.
  class_metrics <- calculate_class_metrics(
    actual = repeat_actual,
    predicted = repeat_predicted,
    group_names = group_names
  )

  class_metrics$Repeat <- repeat_id
  class_metrics$AUC <- as.numeric(
    auc_metrics[paste0("AUC_Class_", 0:(length(group_names) - 1))]
  )

  classwise_repeat_metrics <- rbind(
    classwise_repeat_metrics,
    class_metrics[, c("Repeat", "Class", "Precision", "Recall", "F1", "AUC")]
  )

  repeat_metrics <- rbind(
    repeat_metrics,
    data.frame(
      Repeat = repeat_id,
      t(metrics),
      t(auc_metrics),
      check.names = FALSE
    )
  )
}

# 8. PERFORMANCE SUMMARY

summary_metric_names <- setdiff(names(repeat_metrics), "Repeat")

performance_summary <- data.frame(
  Metric = summary_metric_names,
  Mean = sapply(repeat_metrics[, summary_metric_names, drop = FALSE],
                mean, na.rm = TRUE),
  SD = sapply(repeat_metrics[, summary_metric_names, drop = FALSE],
              sd, na.rm = TRUE)
)

cat("\\nSVM PERFORMANCE SUMMARY INCLUDING AUC:\\n")
print(performance_summary)

# Summary of class-wise metrics across all repeats (mean +/- SD).
classwise_summary <- classwise_repeat_metrics %>%
  group_by(Class) %>%
  summarise(
    Precision_Mean = mean(Precision, na.rm = TRUE),
    Precision_SD = sd(Precision, na.rm = TRUE),
    Recall_Mean = mean(Recall, na.rm = TRUE),
    Recall_SD = sd(Recall, na.rm = TRUE),
    F1_Mean = mean(F1, na.rm = TRUE),
    F1_SD = sd(F1, na.rm = TRUE),
    AUC_Mean = mean(AUC, na.rm = TRUE),
    AUC_SD = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  )

classwise_summary$Class <- factor(classwise_summary$Class, levels = group_names)
classwise_summary <- classwise_summary %>% arrange(Class)

cat("\nCLASS-WISE PERFORMANCE SUMMARY (MEAN +/- SD):\n")
print(classwise_summary)

# 9. POOLED CONFUSION MATRIX

confusion_matrix <- table(
  Actual = factor(
    all_predictions$Actual,
    levels = group_names
  ),
  Predicted = factor(
    all_predictions$Predicted,
    levels = group_names
  )
)

cat("\nPOOLED CONFUSION MATRIX:\n")
print(confusion_matrix)

# Normalize rows to percentages
confusion_normalized <- prop.table(
  confusion_matrix,
  margin = 1
) * 100

confusion_df <- as.data.frame(confusion_normalized)

colnames(confusion_df) <- c(
  "Actual",
  "Predicted",
  "Percentage"
)

# 10. CREATE OUTPUT FOLDER

output_folder <- file.path(
  getwd(),
  "SVM_results"
)

dir.create(
  output_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

# 11. EXPORT RESULTS

write.csv(
  all_predictions,
  file.path(output_folder, "SVM_All_Predictions.csv"),
  row.names = FALSE
)

write.csv(
  repeat_metrics,
  file.path(output_folder, "SVM_Repeat_Metrics.csv"),
  row.names = FALSE
)

write.csv(
  performance_summary,
  file.path(output_folder, "SVM_Performance_Summary.csv"),
  row.names = FALSE
)

write.csv(
  classwise_repeat_metrics,
  file.path(output_folder, "SVM_Classwise_Metrics_All_Repeats.csv"),
  row.names = FALSE
)

write.csv(
  classwise_summary,
  file.path(output_folder, "SVM_Classwise_Performance_Summary.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(confusion_matrix),
  file.path(output_folder, "SVM_Confusion_Matrix.csv")
)

# 12. CONFUSION MATRIX HEATMAP

confusion_plot <- ggplot(
  confusion_df,
  aes(
    x = Predicted,
    y = Actual,
    fill = Percentage
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.1f", Percentage)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue"
  ) +
  labs(
    title = "Linear SVM: Confusion Matrix",
    subtitle = "Row-normalized predictions across repeated cross-validation",
    x = "Predicted group",
    y = "Actual group",
    fill = "Percentage (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_folder, "SVM_Confusion_Matrix.png"),
  plot = confusion_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# 13. PERFORMANCE PLOT

performance_plot <- ggplot(
  performance_summary,
  aes(
    x = Metric,
    y = Mean
  )
) +
  geom_col(fill = "steelblue", width = 0.65) +
  geom_errorbar(
    aes(
      ymin = pmax(0, Mean - SD),
      ymax = pmin(1, Mean + SD)
    ),
    width = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Mean)),
    vjust = -0.6,
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1.1)
  ) +
  labs(
    title = "Linear SVM Classification Performance",
    subtitle = "Mean ± SD across 10 repeated stratified 3-fold CV runs; includes macro-AUC",
    x = NULL,
    y = "Performance score"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(output_folder, "SVM_Performance.png"),
  plot = performance_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Class-wise F1 and AUC plot with mean +/- SD across repeats.
classwise_plot_df <- bind_rows(
  data.frame(
    Class = classwise_summary$Class,
    Metric = "F1-score",
    Mean = classwise_summary$F1_Mean,
    SD = classwise_summary$F1_SD
  ),
  data.frame(
    Class = classwise_summary$Class,
    Metric = "AUC",
    Mean = classwise_summary$AUC_Mean,
    SD = classwise_summary$AUC_SD
  )
)

classwise_plot <- ggplot(
  classwise_plot_df,
  aes(x = Class, y = Mean, fill = Metric)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SD), ymax = pmin(1, Mean + SD)),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "Linear SVM: Class-wise F1-score and AUC",
    subtitle = "Mean +/- SD across 10 repeated stratified 3-fold CV runs",
    x = "Experimental group",
    y = "Performance score",
    fill = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_folder, "SVM_Classwise_F1_AUC.png"),
  plot = classwise_plot,
  width = 11,
  height = 6,
  dpi = 300
)

# 14. FIT FINAL SVM MODEL USING ALL SAMPLES
# This model is for exploratory use, not independent validation.

final_processed <- prepare_data(X, X)

final_svm_model <- svm(
  x = final_processed$train,
  y = group,
  type = "C-classification",
  kernel = "linear",
  cost = 1,
  scale = FALSE
)

saveRDS(
  final_svm_model,
  file.path(output_folder, "Final_Linear_SVM_Model.rds")
)


# 15. PAIRWISE LINEAR SVM PROTEIN WEIGHTS AND TOP 20 PROTEINS
# Exploratory feature-weight analysis using all samples.
# This is separate from cross-validation performance evaluation.
# To match the existing protein-importance outputs, comparisons are:
# Fe3O4_Positive_In_vivo versus each of the other seven groups.

reference_group <- "Fe3O4_Positive_In_vivo"
comparison_groups <- setdiff(group_names, reference_group)

pairwise_weight_list <- list()
top20_list <- list()

for (comparison_group in comparison_groups) {

  cat("\nCalculating protein weights:",
      reference_group, "vs", comparison_group, "\n")

  selected_samples <- group %in% c(reference_group, comparison_group)

  pair_x <- final_processed$train[selected_samples, , drop = FALSE]
  pair_y <- droplevels(factor(
    as.character(group[selected_samples]),
    levels = c(reference_group, comparison_group)
  ))

  pair_model <- svm(
    x = pair_x,
    y = pair_y,
    type = "C-classification",
    kernel = "linear",
    cost = 1,
    scale = FALSE,
    probability = FALSE
  )

  # For a binary linear SVM, the feature-weight vector is:
  # w = sum_i (alpha_i * y_i * support_vector_i).
  # e1071 stores the signed dual coefficients in pair_model$coefs.
  protein_weights <- as.numeric(
    crossprod(pair_model$coefs, pair_model$SV)
  )

  if (length(protein_weights) != ncol(final_processed$train)) {
    stop("Protein-weight vector does not match the feature count.")
  }

  weight_df <- data.frame(
    Protein = colnames(final_processed$train),
    Group_1 = reference_group,
    Group_2 = comparison_group,
    Weight = protein_weights,
    Absolute_Weight = abs(protein_weights),
    stringsAsFactors = FALSE
  )

  weight_df <- weight_df[
    order(weight_df$Absolute_Weight, decreasing = TRUE),
  ]

  weight_df$Rank <- seq_len(nrow(weight_df))

  pairwise_weight_list[[comparison_group]] <- weight_df

  top20_list[[comparison_group]] <- head(weight_df, 20)
}

pairwise_protein_weights <- bind_rows(pairwise_weight_list)
top20_proteins <- bind_rows(top20_list)

# Overall importance across the seven reference-versus-group comparisons.
overall_protein_importance <- pairwise_protein_weights %>%
  group_by(Protein) %>%
  summarise(
    Mean_Absolute_Weight = mean(Absolute_Weight, na.rm = TRUE),
    Maximum_Absolute_Weight = max(Absolute_Weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Absolute_Weight))

# Extract the overall Top 20 proteins ranked by mean absolute SVM weight.
# This produces the four-column table requested for reporting.
overall_top20_proteins <- overall_protein_importance %>%
  slice_head(n = 20) %>%
  mutate(Rank = row_number()) %>%
  transmute(
    Rank = Rank,
    `Protein accession` = Protein,
    `Mean absolute weight` = Mean_Absolute_Weight,
    `Maximum absolute weight` = Maximum_Absolute_Weight
  )

cat("\nOverall Top 20 proteins across all seven pairwise comparisons:\n")
print(overall_top20_proteins, n = 20)

write.csv(
  overall_top20_proteins,
  file.path(output_folder, "SVM_Top20_Overall_Proteins.csv"),
  row.names = FALSE
)

write.csv(
  overall_protein_importance,
  file.path(output_folder, "SVM_Overall_Protein_Importance.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_protein_weights,
  file.path(output_folder, "SVM_Pairwise_Protein_Weights.csv"),
  row.names = FALSE
)

write.csv(
  top20_proteins,
  file.path(output_folder, "SVM_Top20_Proteins_Each_Comparison.csv"),
  row.names = FALSE
)

cat("\nTop 20 proteins for each reference-group comparison:\n")
print(
  top20_proteins %>%
    group_by(Group_2) %>%
    slice_min(order_by = Rank, n = 20, with_ties = FALSE) %>%
    ungroup()
)

cat("\nProtein feature-importance files generated:\n")
cat("1. SVM_Overall_Protein_Importance.csv\n")
cat("2. SVM_Pairwise_Protein_Weights.csv\n")
cat("3. SVM_Top20_Proteins_Each_Comparison.csv\n")
cat("4. SVM_Top20_Overall_Proteins.csv\n")

# 16. FINISH


cat("\n====================================\n")
cat("SVM ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("====================================\n")

cat("\nResults saved in:\n")
cat(output_folder, "\n")

cat("\nFiles generated:\n")
cat("1. SVM_All_Predictions.csv\n")
cat("2. SVM_Repeat_Metrics.csv\n")
cat("3. SVM_Performance_Summary.csv\n")
cat("4. SVM_Confusion_Matrix.csv\n")
cat("5. SVM_Confusion_Matrix.png\n")
cat("6. SVM_Performance.png\n")
cat("7. SVM_Classwise_Metrics_All_Repeats.csv\n")
cat("8. SVM_Classwise_Performance_Summary.csv\n")
cat("9. SVM_Classwise_F1_AUC.png\n")
cat("10. Final_Linear_SVM_Model.rds\n")

# Open output folder
if (.Platform$OS.type == "windows") {
  shell.exec(normalizePath(output_folder))
} else {
  browseURL(output_folder)
}
