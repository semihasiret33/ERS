# ============================================================================
# Script: correction_methods.R
# Purpose: ERS correction methods for cross-country comparisons
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

library(tidyverse)


#' Apply Z-Score Standardization (Person-Mean Centering)
#'
#' Implements Model 1: Z-score standardization by person-mean centering.
#' Each respondent's answers are standardized using their own mean and SD
#' across items, removing person-specific response style effects.
#'
#' @param data Data frame containing response data
#' @param items Character vector of item column names to standardize
#' @param suffix Suffix to add to corrected item names (default = "_zscore")
#' @param keep_original Logical, keep original items? Default is TRUE
#' @return Data frame with z-score standardized items added
#' @references Common method in survey research for ERS correction
#' @examples
#' \dontrun{
#' data_corrected <- apply_zscore_correction(data, items = paste0("Q", 1:10))
#' }
apply_zscore_correction <- function(data, items,
                                   suffix = "_zscore",
                                   keep_original = TRUE) {

  # Validate inputs
  if (!all(items %in% names(data))) {
    stop("Not all specified items found in data")
  }

  # Select item data
  item_data <- data[, items, drop = FALSE]

  # Calculate person mean and SD
  person_mean <- rowMeans(item_data, na.rm = TRUE)
  person_sd <- apply(item_data, 1, sd, na.rm = TRUE)

  # Standardize each item by person-specific mean and SD
  # Z = (X - person_mean) / person_sd
  zscore_data <- sweep(item_data, 1, person_mean, "-")
  zscore_data <- sweep(zscore_data, 1, person_sd, "/")

  # Add suffix to column names
  colnames(zscore_data) <- paste0(items, suffix)

  # Combine with original data
  if (keep_original) {
    result <- cbind(data, zscore_data)
  } else {
    result <- data
    result[, items] <- zscore_data
    colnames(result)[colnames(result) %in% items] <- colnames(zscore_data)
  }

  return(result)
}


#' Calculate Scale Scores with ERS Covariate Control
#'
#' Implements Model 2: Regress scale scores on ERS index and extract residuals.
#' This removes variance attributable to ERS from the scale scores.
#'
#' @param data Data frame containing response data and ERS index
#' @param items Character vector of item column names
#' @param ers_index Name of ERS index variable in data
#' @param method Method for scale score calculation: "mean" or "sum"
#' @param return_residuals Return residuals (TRUE) or predicted values (FALSE)?
#' @return Data frame with ERS-corrected scale score added
#' @examples
#' \dontrun{
#' data <- data %>%
#'   mutate(scale_score_raw = rowMeans(select(., Q1:Q10), na.rm = TRUE)) %>%
#'   add_ers_index(items = paste0("Q", 1:10))
#' data_corrected <- apply_covariate_correction(data, items = paste0("Q", 1:10),
#'                                              ers_index = "ers_index")
#' }
apply_covariate_correction <- function(data, items,
                                      ers_index = "ers_index",
                                      method = "mean",
                                      return_residuals = TRUE) {

  # Calculate raw scale score if not already present
  item_data <- data[, items, drop = FALSE]

  if (method == "mean") {
    scale_score_raw <- rowMeans(item_data, na.rm = TRUE)
  } else if (method == "sum") {
    scale_score_raw <- rowSums(item_data, na.rm = TRUE)
  } else {
    stop("method must be 'mean' or 'sum'")
  }

  # Check if ERS index exists
  if (!ers_index %in% names(data)) {
    stop(paste("ERS index", ers_index, "not found in data"))
  }

  # Regress scale score on ERS index
  # Remove cases with missing values
  complete_cases <- complete.cases(scale_score_raw, data[[ers_index]])

  lm_fit <- lm(scale_score_raw[complete_cases] ~
                data[[ers_index]][complete_cases])

  # Extract residuals or predicted values
  if (return_residuals) {
    # Residuals represent ERS-free scores
    corrected_score <- rep(NA, length(scale_score_raw))
    corrected_score[complete_cases] <- residuals(lm_fit)
    data$scale_score_ers_corrected <- corrected_score
  } else {
    # Predicted values represent ERS effect
    predicted_score <- rep(NA, length(scale_score_raw))
    predicted_score[complete_cases] <- fitted(lm_fit)
    data$scale_score_ers_predicted <- predicted_score
  }

  # Also add raw scale score for comparison
  data$scale_score_raw <- scale_score_raw

  return(data)
}


#' Apply ERS Correction Using Multiple Group IRT
#'
#' Placeholder function for Model 3: ML-MNRM approach.
#' This function would call mirt or TAM to fit multidimensional IRT models
#' with ERS as a latent dimension.
#'
#' @param data Data frame
#' @param items Character vector of item column names
#' @param group_var Name of grouping variable (e.g., country)
#' @param model_type Type of IRT model: "MNRM" or "IRTree"
#' @param ... Additional arguments passed to mirt::multipleGroup()
#' @return List containing fitted model object and extracted parameters
#' @note This is a wrapper function. See separate scripts for full implementation.
apply_mlmnrm_correction <- function(data, items, group_var,
                                   model_type = "MNRM", ...) {

  message("This function requires mirt or TAM package.")
  message("See R/analysis/06_model3_mlmnrm.R for full implementation.")

  # Placeholder - actual implementation in separate script
  stop("Use the full MNRM/IRTree implementation in analysis scripts")
}


#' Compare Country Means Across Correction Methods
#'
#' Calculate and compare country means for raw and corrected scores
#'
#' @param data Data frame with country identifier and scores
#' @param country_var Name of country variable
#' @param score_vars Character vector of score variable names to compare
#' @return Data frame with country means for each score type
#' @examples
#' \dontrun{
#' comparison <- compare_country_means(
#'   data,
#'   country_var = "CNT",
#'   score_vars = c("scale_score_raw", "scale_score_zscore",
#'                  "scale_score_ers_corrected")
#' )
#' }
compare_country_means <- function(data, country_var, score_vars) {

  # Group by country and calculate means
  country_means <- data %>%
    group_by(across(all_of(country_var))) %>%
    summarise(
      across(all_of(score_vars), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      n = n(),
      .groups = "drop"
    )

  return(country_means)
}


#' Calculate Rank-Order Correlation Between Country Rankings
#'
#' Compare how country rankings change across different correction methods
#'
#' @param country_means Data frame with country means (output from
#'   compare_country_means)
#' @param method1 Name of first method's mean column
#' @param method2 Name of second method's mean column
#' @return Spearman rank correlation coefficient
#' @examples
#' \dontrun{
#' cor <- calculate_ranking_correlation(
#'   country_means,
#'   method1 = "mean_scale_score_raw",
#'   method2 = "mean_scale_score_zscore"
#' )
#' }
calculate_ranking_correlation <- function(country_means, method1, method2) {

  # Extract the two columns
  if (!method1 %in% names(country_means) | !method2 %in% names(country_means)) {
    stop("Specified method columns not found in data")
  }

  # Calculate Spearman correlation
  cor_result <- cor(country_means[[method1]],
                   country_means[[method2]],
                   method = "spearman",
                   use = "complete.obs")

  return(cor_result)
}


#' Create Long-Format Comparison Table
#'
#' Reshape country means for easier comparison and visualization
#'
#' @param country_means Data frame with country means (wide format)
#' @param country_var Name of country variable
#' @param value_vars Character vector of mean columns to reshape
#' @return Data frame in long format
create_comparison_table <- function(country_means, country_var, value_vars) {

  # Reshape to long format
  long_data <- country_means %>%
    select(all_of(c(country_var, value_vars))) %>%
    pivot_longer(
      cols = all_of(value_vars),
      names_to = "method",
      values_to = "country_mean"
    ) %>%
    # Clean method names (remove "mean_" prefix if present)
    mutate(method = str_remove(method, "^mean_"))

  return(long_data)
}
