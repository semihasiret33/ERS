# ============================================================================
# Script: 04_ers_classification_accuracy.R
# Purpose: Evaluate ERS detection accuracy for each correction method
# Author: ERS Research Project
# Date: 2025-11-30
# ============================================================================
#
# This script evaluates how well each method identifies individuals with
# high vs. low Extreme Response Style (ERS).
#
# Metrics calculated:
# - Sensitivity: Of truly high-ERS individuals, % correctly identified
# - Specificity: Of truly low-ERS individuals, % correctly identified
# - Overall Accuracy: % of all individuals correctly classified
# - Positive Predictive Value (Precision)
# - F1 Score
#
# Classification threshold: Median ERS value in the sample
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(caret)  # For confusion matrix metrics

cat("\n")
cat("=======================================================\n")
cat("  ERS CLASSIFICATION ACCURACY ANALYSIS                \n")
cat("=======================================================\n")
cat("\n")

# 2. Load Simulated Data with Correction Results ----

cat("Loading simulation data...\n")

# Load original data with true ERS values
all_datasets <- readRDS("data/simulated/all_conditions_full_factorial.rds")

cat("  Data loaded successfully\n")
cat("  Total conditions:", length(all_datasets), "\n")
cat("\n")


# 3. Helper Functions ----

#' Classify ERS as High or Low Based on Threshold
#'
#' @param ers_values Numeric vector of ERS values
#' @param threshold Classification threshold (default = median)
#' @return Binary vector (1 = high ERS, 0 = low ERS)
classify_ers <- function(ers_values, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- median(ers_values, na.rm = TRUE)
  }

  classification <- ifelse(ers_values > threshold, 1, 0)
  return(classification)
}


#' Calculate Classification Metrics
#'
#' @param true_class True binary classification (1/0)
#' @param pred_class Predicted binary classification (1/0)
#' @return Data frame with accuracy metrics
calculate_classification_metrics <- function(true_class, pred_class) {

  # Remove NAs
  valid_cases <- complete.cases(true_class, pred_class)
  true_class <- true_class[valid_cases]
  pred_class <- pred_class[valid_cases]

  # Create confusion matrix
  # TRUE POSITIVES: Truly high ERS, predicted high ERS
  tp <- sum(true_class == 1 & pred_class == 1)

  # TRUE NEGATIVES: Truly low ERS, predicted low ERS
  tn <- sum(true_class == 0 & pred_class == 0)

  # FALSE POSITIVES: Truly low ERS, predicted high ERS
  fp <- sum(true_class == 0 & pred_class == 1)

  # FALSE NEGATIVES: Truly high ERS, predicted low ERS
  fn <- sum(true_class == 1 & pred_class == 0)

  # Calculate metrics
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA  # True Positive Rate
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA  # True Negative Rate
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA  # Positive Predictive Value

  # F1 Score (harmonic mean of precision and recall)
  f1_score <- if (!is.na(precision) && !is.na(sensitivity)) {
    2 * (precision * sensitivity) / (precision + sensitivity)
  } else {
    NA
  }

  # Return metrics
  metrics <- data.frame(
    n_total = length(true_class),
    n_high_ers_true = sum(true_class == 1),
    n_low_ers_true = sum(true_class == 0),
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision,
    f1_score = f1_score
  )

  return(metrics)
}


#' Estimate ERS from Corrected Scores
#'
#' For each correction method, we need to estimate ERS.
#' Since methods remove ERS, we look at the residuals or compare to raw.
#'
#' @param raw_scores Original scale scores
#' @param corrected_scores Corrected scale scores
#' @return Estimated ERS index (higher = more ERS)
estimate_ers_from_correction <- function(raw_scores, corrected_scores) {

  # ERS estimate: Difference between raw and corrected scores
  # Larger difference suggests more ERS was removed
  ers_estimate <- abs(raw_scores - corrected_scores)

  return(ers_estimate)
}


# 4. Analyze ERS Classification for All Conditions ----

cat("=== Analyzing ERS Classification Accuracy ===\n")
cat("\n")

all_classification_results <- list()

# Select a subset of conditions for testing (can expand to all 72 later)
# For now, test with first 12 conditions to ensure it works
test_conditions <- names(all_datasets)[1:min(12, length(all_datasets))]

for (condition_name in test_conditions) {

  cat("\n--- Processing:", condition_name, "---\n")

  replications <- all_datasets[[condition_name]]
  n_reps <- length(replications)

  # Process a subset of replications (e.g., first 10)
  n_reps_to_process <- min(10, n_reps)

  condition_metrics <- list()

  for (rep in 1:n_reps_to_process) {

    if (rep %% 5 == 0) {
      cat("  Replication", rep, "/", n_reps_to_process, "\n")
    }

    # Get data
    sim_data <- replications[[rep]]$data
    true_ers <- sim_data$true_ers

    # Classify true ERS (high vs low based on median)
    true_ers_class <- classify_ers(true_ers)

    # Get item columns
    item_cols <- grep("^Item", names(sim_data), value = TRUE)

    # Calculate raw scores
    raw_scores <- rowMeans(sim_data[, item_cols], na.rm = TRUE)

    # Method 1: Z-score
    # Apply z-score correction
    source("R/functions/correction_methods.R")
    source("R/functions/ers_indices.R")

    data_zscore <- apply_zscore_correction(sim_data, items = item_cols, suffix = "_z")
    zscore_cols <- paste0(item_cols, "_z")
    zscore_scores <- rowMeans(data_zscore[, zscore_cols], na.rm = TRUE)

    # Estimate ERS from z-score (based on amount of correction)
    ers_est_zscore <- estimate_ers_from_correction(raw_scores, zscore_scores)
    pred_zscore_class <- classify_ers(ers_est_zscore)

    # Calculate metrics for z-score
    metrics_zscore <- calculate_classification_metrics(true_ers_class, pred_zscore_class)
    metrics_zscore$method <- "zscore"
    metrics_zscore$condition <- condition_name
    metrics_zscore$replication <- rep


    # Method 2: Covariate (using actual ERS index)
    data_ers <- add_ers_index(sim_data, items = item_cols, index_name = "ers_index")
    pred_covariate_class <- classify_ers(data_ers$ers_index)

    # Calculate metrics for covariate (ERS index)
    metrics_covariate <- calculate_classification_metrics(true_ers_class, pred_covariate_class)
    metrics_covariate$method <- "covariate"
    metrics_covariate$condition <- condition_name
    metrics_covariate$replication <- rep


    # Store results
    condition_metrics[[paste0("rep", rep, "_zscore")]] <- metrics_zscore
    condition_metrics[[paste0("rep", rep, "_covariate")]] <- metrics_covariate
  }

  # Combine metrics for this condition
  all_classification_results[[condition_name]] <- bind_rows(condition_metrics)
}


# Combine all results
cat("\n\nCombining all classification results...\n")
classification_results <- bind_rows(all_classification_results)

cat("  Total rows:", nrow(classification_results), "\n")
cat("\n")


# 5. Aggregate Results by Method and Condition ----

cat("=== Aggregating Classification Metrics ===\n")
cat("\n")

aggregated_metrics <- classification_results %>%
  group_by(condition, method) %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_f1_score = mean(f1_score, na.rm = TRUE),
    n_replications = n(),
    .groups = "drop"
  )

cat("Sample of aggregated metrics:\n")
print(head(aggregated_metrics, 10))
cat("\n")


# 6. Save Results ----

cat("=== Saving Results ===\n")
cat("\n")

# Save detailed results
write_csv(classification_results,
         file = "output/tables/ers_classification_detailed.csv")
cat("  ✓ Detailed results: output/tables/ers_classification_detailed.csv\n")

# Save aggregated results
write_csv(aggregated_metrics,
         file = "output/tables/ers_classification_summary.csv")
cat("  ✓ Summary results: output/tables/ers_classification_summary.csv\n")


# 7. Summary Report ----

cat("\n")
cat("=======================================================\n")
cat("  ERS CLASSIFICATION ACCURACY SUMMARY                 \n")
cat("=======================================================\n")
cat("\n")

cat("Overall Performance by Method:\n")
overall_performance <- aggregated_metrics %>%
  group_by(method) %>%
  summarise(
    avg_sensitivity = mean(mean_sensitivity, na.rm = TRUE),
    avg_specificity = mean(mean_specificity, na.rm = TRUE),
    avg_accuracy = mean(mean_accuracy, na.rm = TRUE),
    avg_f1_score = mean(mean_f1_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_accuracy))

print(overall_performance)

cat("\n")
cat("Interpretation:\n")
cat("  - Sensitivity: Yüksek ERS'ye sahip bireylerin doğru tespit oranı\n")
cat("  - Specificity: Düşük ERS'ye sahip bireylerin doğru tespit oranı\n")
cat("  - Accuracy: Genel doğruluk (tüm bireylerin doğru sınıflandırılma oranı)\n")
cat("  - F1 Score: Hassasiyet ve duyarlılığın harmonik ortalaması\n")
cat("\n")

best_method <- overall_performance$method[1]
cat("En iyi performans gösteren yöntem:", best_method, "\n")
cat("  Accuracy:", round(overall_performance$avg_accuracy[1] * 100, 2), "%\n")
cat("\n")

cat("=======================================================\n")
cat("Script completed:", as.character(Sys.time()), "\n")
cat("=======================================================\n")
