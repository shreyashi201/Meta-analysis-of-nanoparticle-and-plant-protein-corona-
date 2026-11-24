
# ============================================================
# RANDOM FOREST CLASSIFICATION OF PROTEOMICS DATA
# 8 GROUPS | 3 BIOLOGICAL REPLICATES | 24 SAMPLES
# ============================================================


# STEP 1: INSTALL REQUIRED PACKAGES

if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest",
                   repos = "https://cloud.r-project.org")
}

library(randomForest)

if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC", repos = "https://cloud.r-project.org")
}
library(pROC)


# STEP 2: IMPORT DATASET

dat <- read.csv(
  file.choose(),
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# STEP 3: EXTRACT PROTEIN IDs AND LABELS

protein_ids <- as.character(dat[-1, 1])

labels <- as.integer(
  as.character(unlist(dat[1, -1]))
)


# STEP 4: EXTRACT PROTEIN ABUNDANCE DATA

protein_data <- dat[-1, -1]

protein_data[] <- lapply(
  protein_data,
  function(x) as.numeric(as.character(x))
)

protein_matrix <- as.matrix(protein_data)


# STEP 5: TRANSPOSE DATA
# Rows = Samples
# Columns = Proteins

X <- as.data.frame(t(protein_matrix))

colnames(X) <- make.unique(protein_ids)

rownames(X) <- colnames(dat)[-1]


# STEP 6: DEFINE EXPERIMENTAL GROUPS

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

y <- factor(
  labels,
  levels = 0:7,
  labels = group_names
)


# STEP 7: CHECK DATA

cat("Number of samples:", nrow(X), "\n")
cat("Number of proteins:", ncol(X), "\n")

print(table(y))

if (anyNA(y)) {
  stop("Invalid or missing group labels detected.")
}

if (anyNA(X)) {
  stop("Missing values detected. Handle missing values before modelling.")
}

if (any(!is.finite(as.matrix(X)))) {
  stop("Non-finite abundance values detected.")
}


# STEP 8: REMOVE NON-INFORMATIVE PROTEINS

keep <- apply(
  X,
  2,
  function(z) length(unique(z)) > 1
)

X <- X[, keep, drop = FALSE]

cat("Proteins retained:", ncol(X), "\n")


# STEP 9: CHECK GROUP REPLICATES

if (any(table(y) < 3)) {
  stop("At least three samples per group are required for this 3-fold CV.")
}


# STEP 10: CREATE STRATIFIED 3-FOLD CROSS-VALIDATION

set.seed(42)

number_repeats <- 10
number_folds <- 3

cv_results <- data.frame()

all_predictions <- data.frame()

# One-vs-rest AUC for all eight classes; macro-AUC is their unweighted mean.
calculate_ovr_auc <- function(actual, probability_matrix, class_levels) {
  auc_values <- setNames(rep(NA_real_, length(class_levels)),
                         paste0("AUC_Class_", seq_along(class_levels)-1))

  for (i in seq_along(class_levels)) {
    target_class <- class_levels[i]
    binary_actual <- as.integer(as.character(actual) == target_class)

    if (length(unique(binary_actual)) == 2 &&
        target_class %in% colnames(probability_matrix)) {
      roc_obj <- tryCatch(
        pROC::roc(response = binary_actual,
                  predictor = probability_matrix[, target_class],
                  levels = c(0, 1),
                  direction = "<",
                  quiet = TRUE),
        error = function(e) NULL
      )
      if (!is.null(roc_obj)) {
        auc_values[i] <- as.numeric(pROC::auc(roc_obj))
      }
    }
  }

  c(Macro_AUC = mean(auc_values, na.rm = TRUE), auc_values)
}


# STEP 11: REPEATED CROSS-VALIDATION

for (r in 1:number_repeats) {
  
  set.seed(100 + r)
  
  # Randomly assign samples to three folds within each group
  fold_assignment <- unlist(
    lapply(levels(y), function(g) {
      
      group_indices <- which(y == g)
      
      sample(rep(1:number_folds,
                 length.out = length(group_indices)))
    })
  )
  
  names(fold_assignment) <- NULL
  
  # Reconstruct assignments in original sample order
  fold_id <- integer(length(y))
  
  position <- 1
  
  for (g in levels(y)) {
    
    group_indices <- which(y == g)
    
    fold_id[group_indices] <-
      fold_assignment[position:(position + length(group_indices) - 1)]
    
    position <- position + length(group_indices)
  }
  
  
  for (f in 1:number_folds) {
    
    train_indices <- which(fold_id != f)
    test_indices <- which(fold_id == f)
    
    X_train <- X[train_indices, , drop = FALSE]
    X_test <- X[test_indices, , drop = FALSE]
    
    y_train <- droplevels(y[train_indices])
    y_test <- y[test_indices]
    
    
    # Random Forest model
    set.seed(1000 + r * 10 + f)
    
    rf_model <- randomForest(
      x = X_train,
      y = y_train,
      ntree = 1000,
      mtry = floor(sqrt(ncol(X_train))),
      importance = TRUE
    )
    
    
    # Predicted classes and class probabilities
    predictions <- predict(
      rf_model,
      newdata = X_test,
      type = "response"
    )

    probabilities <- predict(
      rf_model,
      newdata = X_test,
      type = "prob"
    )
    
    
    # Confusion matrix
    cm <- table(
      factor(y_test, levels = levels(y)),
      factor(predictions, levels = levels(y))
    )
    
    
    # Metrics
    TP <- diag(cm)
    
    precision <- TP / colSums(cm)
    recall <- TP / rowSums(cm)
    
    F1 <- 2 * precision * recall /
      (precision + recall)
    
    accuracy <- sum(TP) / sum(cm)
    
    macro_precision <- mean(precision, na.rm = TRUE)
    macro_recall <- mean(recall, na.rm = TRUE)
    macro_F1 <- mean(F1, na.rm = TRUE)
    
    
    # Store results
    cv_results <- rbind(
      cv_results,
      data.frame(
        Repeat = r,
        Fold = f,
        Accuracy = accuracy,
        Precision = macro_precision,
        Recall = macro_recall,
        F1 = macro_F1
      )
    )
    
    
    # Store individual predictions
    all_predictions <- rbind(
      all_predictions,
      data.frame(
        Sample = rownames(X)[test_indices],
        Actual = as.character(y_test),
        Predicted = as.character(predictions),
        Repeat = r,
        Fold = f,
        as.data.frame(probabilities, check.names = FALSE)
      )
    )
    
  }
}


# STEP 11B: REPEAT-WISE OUT-OF-FOLD METRICS AND MULTICLASS AUC
# Pool the held-out predictions from all three folds within each repeat.

repeat_results <- do.call(rbind, lapply(1:number_repeats, function(r) {
  d <- all_predictions[all_predictions$Repeat == r, , drop = FALSE]

  repeat_cm <- table(
    factor(d$Actual, levels = group_names),
    factor(d$Predicted, levels = group_names)
  )
  tp <- diag(repeat_cm)
  fp <- colSums(repeat_cm) - tp
  fn <- rowSums(repeat_cm) - tp
  class_precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  class_recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  class_f1 <- ifelse(2 * tp + fp + fn == 0, 0,
                     2 * tp / (2 * tp + fp + fn))
  metric_values <- c(
    Accuracy = sum(tp) / sum(repeat_cm),
    Precision = mean(class_precision),
    Recall = mean(class_recall),
    F1 = mean(class_f1)
  )
  probability_columns <- group_names
  probability_matrix <- as.matrix(d[, probability_columns, drop = FALSE])
  storage.mode(probability_matrix) <- "numeric"

  auc_values <- calculate_ovr_auc(
    actual = d$Actual,
    probability_matrix = probability_matrix,
    class_levels = group_names
  )

  data.frame(
    Repeat = r,
    t(metric_values),
    t(auc_values),
    check.names = FALSE
  )
}))

write.csv(repeat_results, "RF_Repeat_Metrics_with_AUC.csv",
          row.names = FALSE)

repeat_metric_names <- setdiff(names(repeat_results), "Repeat")
repeat_performance_summary <- data.frame(
  Metric = repeat_metric_names,
  Mean = sapply(repeat_results[, repeat_metric_names, drop = FALSE],
                mean, na.rm = TRUE),
  SD = sapply(repeat_results[, repeat_metric_names, drop = FALSE],
              sd, na.rm = TRUE)
)

write.csv(repeat_performance_summary,
          "RF_Performance_Summary_with_AUC.csv",
          row.names = FALSE)

cat("\\n========== RANDOM FOREST REPEAT-WISE PERFORMANCE WITH AUC ==========\\n")
print(repeat_performance_summary)


# STEP 12: SUMMARIZE MODEL PERFORMANCE

cat("\n========== RANDOM FOREST PERFORMANCE ==========\n")

performance_summary <- data.frame(
  Metric = c("Accuracy", "Precision", "Recall", "F1-score"),
  
  Mean = c(
    mean(cv_results$Accuracy),
    mean(cv_results$Precision),
    mean(cv_results$Recall),
    mean(cv_results$F1)
  ),
  
  SD = c(
    sd(cv_results$Accuracy),
    sd(cv_results$Precision),
    sd(cv_results$Recall),
    sd(cv_results$F1)
  )
)

print(performance_summary)


# STEP 13: COMBINED OUT-OF-FOLD CONFUSION MATRIX

combined_cm <- table(
  Actual = factor(
    all_predictions$Actual,
    levels = levels(y)
  ),
  Predicted = factor(
    all_predictions$Predicted,
    levels = levels(y)
  )
)

print(combined_cm)


# STEP 14: TRAIN FINAL MODEL ON ALL SAMPLES
# This model is for exploratory protein importance,
# not independent performance estimation.

set.seed(42)

final_rf <- randomForest(
  x = X,
  y = y,
  ntree = 1000,
  mtry = floor(sqrt(ncol(X))),
  importance = TRUE
)


# STEP 15: EXTRACT PROTEIN IMPORTANCE

protein_importance <- importance(
  final_rf,
  type = 2
)

importance_table <- data.frame(
  Protein = rownames(protein_importance),
  Importance = protein_importance[, 1]
)

importance_table <- importance_table[
  order(-importance_table$Importance),
]


# Display top 20 proteins

cat("\n========== TOP 20 IMPORTANT PROTEINS ==========\n")

print(head(importance_table, 20))


# STEP 16: SAVE RESULTS

write.csv(
  performance_summary,
  "RF_Performance_Summary.csv",
  row.names = FALSE
)

write.csv(
  cv_results,
  "RF_CV_Results.csv",
  row.names = FALSE
)

write.csv(
  combined_cm,
  "RF_Confusion_Matrix.csv"
)

write.csv(
  all_predictions,
  "RF_CV_Predictions.csv",
  row.names = FALSE
)

write.csv(
  importance_table,
  "RF_Protein_Importance.csv",
  row.names = FALSE
)


cat("\nAnalysis completed successfully!\n")



# ============================================================
# COMPLETE RANDOM FOREST RESULTS AND VISUALIZATION
# ============================================================

# 1. PERFORMANCE SUMMARY

cat("\n========== MODEL PERFORMANCE ==========\n")
print(performance_summary)


# 2. CONFUSION MATRIX

cat("\n========== CONFUSION MATRIX ==========\n")
print(combined_cm)


# 3. TOP 30 IMPORTANT PROTEINS

cat("\n========== TOP 30 PROTEINS ==========\n")

top30 <- head(importance_table, 30)
print(top30)


# 4. PLOT TOP 20 PROTEINS

top20 <- head(importance_table, 20)

png("RF_Top20_Protein_Importance.png",
    width = 1400, height = 1000, res = 150)

par(mar = c(5, 12, 4, 2))

barplot(
  rev(top20$Importance),
  names.arg = rev(top20$Protein),
  horiz = TRUE,
  las = 1,
  col = "steelblue",
  main = "Top 20 Important Proteins",
  xlab = "Mean Decrease in Gini"
)

dev.off()


# 5. CROSS-VALIDATION PERFORMANCE PLOT

png("RF_Model_Performance.png",
    width = 1200, height = 900, res = 150)

boxplot(
  cv_results[, c("Accuracy", "Precision", "Recall", "F1")],
  col = "lightblue",
  ylab = "Performance",
  main = "Random Forest Cross-validation Performance",
  ylim = c(0, 1)
)

dev.off()

# 5B. REPEAT-WISE PERFORMANCE INCLUDING MACRO-AUC
png("RF_Repeat_Performance_with_AUC.png",
    width = 1400, height = 900, res = 150)

boxplot(
  repeat_results[, c("Accuracy", "Precision", "Recall", "F1", "Macro_AUC")],
  ylab = "Performance",
  main = "Random Forest Repeat-wise Performance",
  ylim = c(0, 1),
  las = 2
)

dev.off()


# 6. CONFUSION MATRIX HEATMAP

png("RF_Confusion_Matrix.png",
    width = 1200, height = 1000, res = 150)

heatmap(
  as.matrix(combined_cm),
  Rowv = NA,
  Colv = NA,
  scale = "none",
  margins = c(10, 10),
  main = "Random Forest Confusion Matrix"
)

dev.off()


# 7. SAVE RESULTS

write.csv(
  performance_summary,
  "RF_Performance_Summary.csv",
  row.names = FALSE
)

write.csv(
  top30,
  "RF_Top30_Proteins.csv",
  row.names = FALSE
)

write.csv(
  combined_cm,
  "RF_Confusion_Matrix.csv"
)

write.csv(
  cv_results,
  "RF_CV_Results.csv",
  row.names = FALSE
)

write.csv(
  all_predictions,
  "RF_CV_Predictions.csv",
  row.names = FALSE
)


# 8. FINAL MESSAGE

cat("\n====================================\n")
cat("RANDOM FOREST ANALYSIS COMPLETED\n")
cat("====================================\n")

cat("\nPerformance summary with macro-AUC, repeat-wise metrics,")
cat("\nconfusion matrix, protein importance and plots saved.")

cat("\n\nFiles are saved in:\n")
print(getwd())

# ============================================================
# CORRECTED CROSS-VALIDATION PERFORMANCE
# ============================================================

groups <- levels(y)

calculate_metrics <- function(actual, predicted) {
  
  cm <- table(
    factor(actual, levels = groups),
    factor(predicted, levels = groups)
  )
  
  TP <- diag(cm)
  FP <- colSums(cm) - TP
  FN <- rowSums(cm) - TP
  
  precision <- ifelse(TP + FP == 0, 0,
                      TP / (TP + FP))
  
  recall <- ifelse(TP + FN == 0, 0,
                   TP / (TP + FN))
  
  f1 <- ifelse(
    2 * TP + FP + FN == 0,
    0,
    2 * TP / (2 * TP + FP + FN)
  )
  
  c(
    Accuracy = sum(TP) / sum(cm),
    Precision = mean(precision),
    Recall = mean(recall),
    F1 = mean(f1)
  )
}

# Calculate metrics separately for every fold
corrected_results <- do.call(
  rbind,
  lapply(
    split(all_predictions,
          interaction(all_predictions$Repeat,
                      all_predictions$Fold,
                      drop = TRUE)),
    function(d) {
      calculate_metrics(d$Actual, d$Predicted)
    }
  )
)

# Mean and SD across folds
corrected_summary <- data.frame(
  Metric = colnames(corrected_results),
  Mean = colMeans(corrected_results),
  SD = apply(corrected_results, 2, sd)
)

print(corrected_summary)

write.csv(
  corrected_summary,
  "RF_Corrected_Performance.csv",
  row.names = FALSE
)

