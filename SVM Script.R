
############################################################
# MULTICLASS LINEAR SVM ANALYSIS
# Dataset: 8 experimental groups, 3 replicates/group
# Method: Repeated stratified 3-fold cross-validation
# Repeats: 10
############################################################

# 1. INSTALL AND LOAD PACKAGES
packages <- c("e1071", "ggplot2", "dplyr", "tidyr")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)

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
  
  return(c(
    Accuracy = accuracy,
    Macro_Precision = mean(precision),
    Macro_Recall = mean(recall),
    Macro_F1 = mean(f1)
  ))
}

# 7. REPEATED STRATIFIED 3-FOLD CROSS-VALIDATION

n_folds <- 3
n_repeats <- 10

all_predictions <- data.frame()
repeat_metrics <- data.frame()

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
      decision.values = FALSE
    )
    
    predictions <- predict(
      svm_model,
      newdata = processed$test
    )
    
    predictions <- factor(
      predictions,
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
  
  # Calculate metrics for this complete repeat
  metrics <- calculate_metrics(
    repeat_actual,
    repeat_predicted,
    group_names
  )
  
  repeat_metrics <- rbind(
    repeat_metrics,
    data.frame(
      Repeat = repeat_id,
      t(metrics)
    )
  )
}

# 8. PERFORMANCE SUMMARY

performance_summary <- data.frame(
  Metric = c(
    "Accuracy",
    "Macro Precision",
    "Macro Recall",
    "Macro F1"
  ),
  Mean = c(
    mean(repeat_metrics$Accuracy),
    mean(repeat_metrics$Macro_Precision),
    mean(repeat_metrics$Macro_Recall),
    mean(repeat_metrics$Macro_F1)
  ),
  SD = c(
    sd(repeat_metrics$Accuracy),
    sd(repeat_metrics$Macro_Precision),
    sd(repeat_metrics$Macro_Recall),
    sd(repeat_metrics$Macro_F1)
  )
)

cat("\nSVM PERFORMANCE SUMMARY:\n")
print(performance_summary)

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
    subtitle = "Mean ± SD across 10 repeated stratified 3-fold CV runs",
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

# 15. FINISH

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
cat("7. Final_Linear_SVM_Model.rds\n")

# Open output folder
if (.Platform$OS.type == "windows") {
  shell.exec(normalizePath(output_folder))
} else {
  browseURL(output_folder)
}

