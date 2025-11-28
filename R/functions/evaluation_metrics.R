# ============================================================================
# Script: evaluation_metrics.R
# Purpose: Functions to calculate evaluation metrics (Bias, RMSE, etc.)
#          for comparing ERS correction methods
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

#' Calculate Bias
#'
#' Computes the bias between estimated and true values.
#' Bias = mean(estimate - true_value)
#'
#' @param estimates Numeric vector of estimated values
#' @param true_values Numeric vector of true values
#' @param na_rm Logical, remove NA values? Default is TRUE
#' @return Numeric, the bias value
#' @examples
#' \dontrun{
#' bias <- calculate_bias(estimates = estimated_means, true_values = true_means)
#' }
calculate_bias <- function(estimates, true_values, na_rm = TRUE) {

  if (length(estimates) != length(true_values)) {
    stop("estimates and true_values must have the same length")
  }

  bias <- mean(estimates - true_values, na.rm = na_rm)

  return(bias)
}


#' Calculate Root Mean Square Error (RMSE)
#'
#' Computes RMSE between estimated and true values.
#' RMSE = sqrt(mean((estimate - true_value)^2))
#'
#' @param estimates Numeric vector of estimated values
#' @param true_values Numeric vector of true values
#' @param na_rm Logical, remove NA values? Default is TRUE
#' @return Numeric, the RMSE value
#' @examples
#' \dontrun{
#' rmse <- calculate_rmse(estimates = estimated_means, true_values = true_means)
#' }
calculate_rmse <- function(estimates, true_values, na_rm = TRUE) {

  if (length(estimates) != length(true_values)) {
    stop("estimates and true_values must have the same length")
  }

  squared_errors <- (estimates - true_values)^2
  rmse <- sqrt(mean(squared_errors, na.rm = na_rm))

  return(rmse)
}


#' Calculate Mean Absolute Error (MAE)
#'
#' Computes MAE between estimated and true values.
#' MAE = mean(|estimate - true_value|)
#'
#' @param estimates Numeric vector of estimated values
#' @param true_values Numeric vector of true values
#' @param na_rm Logical, remove NA values? Default is TRUE
#' @return Numeric, the MAE value
calculate_mae <- function(estimates, true_values, na_rm = TRUE) {

  if (length(estimates) != length(true_values)) {
    stop("estimates and true_values must have the same length")
  }

  absolute_errors <- abs(estimates - true_values)
  mae <- mean(absolute_errors, na.rm = na_rm)

  return(mae)
}


#' Calculate Relative Bias
#'
#' Computes relative bias as a percentage.
#' Relative Bias = (Bias / mean(true_values)) * 100
#'
#' @param estimates Numeric vector of estimated values
#' @param true_values Numeric vector of true values
#' @param na_rm Logical, remove NA values? Default is TRUE
#' @return Numeric, the relative bias in percentage
calculate_relative_bias <- function(estimates, true_values, na_rm = TRUE) {

  bias <- calculate_bias(estimates, true_values, na_rm)
  mean_true <- mean(true_values, na.rm = na_rm)

  if (mean_true == 0) {
    warning("Mean of true values is zero, relative bias is undefined")
    return(NA)
  }

  relative_bias <- (bias / mean_true) * 100

  return(relative_bias)
}


#' Evaluate Correction Method Performance
#'
#' Comprehensive evaluation of a correction method's performance across
#' multiple metrics
#'
#' @param estimates Numeric vector of estimated country means
#' @param true_values Numeric vector of true country means
#' @param method_name Character, name of the correction method
#' @return Data frame with evaluation metrics
#' @examples
#' \dontrun{
#' results <- evaluate_method(
#'   estimates = country_means_model1,
#'   true_values = true_country_means,
#'   method_name = "Z-score"
#' )
#' }
evaluate_method <- function(estimates, true_values, method_name = "Method") {

  results <- data.frame(
    method = method_name,
    bias = calculate_bias(estimates, true_values),
    rmse = calculate_rmse(estimates, true_values),
    mae = calculate_mae(estimates, true_values),
    relative_bias = calculate_relative_bias(estimates, true_values),
    correlation = cor(estimates, true_values, use = "complete.obs"),
    n = sum(complete.cases(estimates, true_values)),
    stringsAsFactors = FALSE
  )

  return(results)
}


#' Compare Multiple Correction Methods
#'
#' Evaluate and compare multiple correction methods simultaneously
#'
#' @param estimates_list Named list of estimate vectors
#' @param true_values Numeric vector of true values (same for all methods)
#' @return Data frame with comparison of all methods
#' @examples
#' \dontrun{
#' comparison <- compare_methods(
#'   estimates_list = list(
#'     "Raw" = country_means_raw,
#'     "Z-score" = country_means_zscore,
#'     "Covariate" = country_means_covariate,
#'     "ML-MNRM" = country_means_mlmnrm
#'   ),
#'   true_values = true_country_means
#' )
#' }
compare_methods <- function(estimates_list, true_values) {

  # Check that estimates_list is a named list
  if (is.null(names(estimates_list))) {
    stop("estimates_list must be a named list")
  }

  # Evaluate each method
  results_list <- lapply(names(estimates_list), function(method_name) {
    evaluate_method(
      estimates = estimates_list[[method_name]],
      true_values = true_values,
      method_name = method_name
    )
  })

  # Combine results
  results <- do.call(rbind, results_list)

  # Order by RMSE (best performing first)
  results <- results[order(results$rmse), ]

  return(results)
}


#' Calculate Coverage Probability (for simulation studies)
#'
#' Calculates the proportion of confidence intervals that contain the true value
#'
#' @param estimates Matrix where each column is a replication
#' @param true_value Numeric, the true parameter value
#' @param se_estimates Matrix of standard errors (same dimensions as estimates)
#' @param alpha Significance level (default = 0.05 for 95% CI)
#' @return Numeric, coverage probability
calculate_coverage <- function(estimates, true_value, se_estimates, alpha = 0.05) {

  # Calculate critical value
  z_crit <- qnorm(1 - alpha/2)

  # Calculate confidence intervals
  ci_lower <- estimates - z_crit * se_estimates
  ci_upper <- estimates + z_crit * se_estimates

  # Check if true value is within CI
  covered <- (ci_lower <= true_value) & (true_value <= ci_upper)

  # Calculate coverage probability
  coverage <- mean(covered, na.rm = TRUE)

  return(coverage)
}


#' Aggregate Simulation Results
#'
#' Summarize results across multiple simulation replications
#'
#' @param sim_results Data frame with columns: replication, method, estimate,
#'   true_value
#' @param group_vars Character vector of grouping variables (e.g., c("method",
#'   "condition"))
#' @return Data frame with aggregated bias and RMSE by group
#' @examples
#' \dontrun{
#' summary <- aggregate_simulation_results(
#'   sim_results,
#'   group_vars = c("method", "ers_level")
#' )
#' }
aggregate_simulation_results <- function(sim_results, group_vars = "method") {

  library(tidyverse)

  summary <- sim_results %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_estimate = mean(estimate, na.rm = TRUE),
      mean_true_value = mean(true_value, na.rm = TRUE),
      bias = mean(estimate - true_value, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_value)^2, na.rm = TRUE)),
      mae = mean(abs(estimate - true_value), na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      n_replications = n(),
      .groups = "drop"
    )

  return(summary)
}


#' Calculate Percent Change in Country Means
#'
#' Calculate how much country means change after ERS correction
#'
#' @param means_raw Numeric vector of raw country means
#' @param means_corrected Numeric vector of corrected country means
#' @return Numeric vector of percent change
calculate_percent_change <- function(means_raw, means_corrected) {

  if (length(means_raw) != length(means_corrected)) {
    stop("Input vectors must have same length")
  }

  percent_change <- ((means_corrected - means_raw) / means_raw) * 100

  return(percent_change)
}


#' Calculate Change in Country Rankings
#'
#' Quantify how much country rankings change between methods
#'
#' @param means_method1 Numeric vector of country means from method 1
#' @param means_method2 Numeric vector of country means from method 2
#' @param country_names Character vector of country names (optional)
#' @return Data frame showing ranking changes
calculate_ranking_change <- function(means_method1, means_method2,
                                    country_names = NULL) {

  if (is.null(country_names)) {
    country_names <- paste0("Country_", 1:length(means_method1))
  }

  # Calculate rankings (1 = highest mean)
  rank_method1 <- rank(-means_method1, ties.method = "average")
  rank_method2 <- rank(-means_method2, ties.method = "average")

  # Calculate rank change
  rank_change <- rank_method2 - rank_method1

  # Create result data frame
  result <- data.frame(
    country = country_names,
    mean_method1 = means_method1,
    rank_method1 = rank_method1,
    mean_method2 = means_method2,
    rank_method2 = rank_method2,
    rank_change = rank_change,
    abs_rank_change = abs(rank_change),
    stringsAsFactors = FALSE
  )

  # Order by absolute rank change (largest changes first)
  result <- result[order(-result$abs_rank_change), ]

  return(result)
}
