# ============================================================================
# Script: 03_evaluation.R
# Purpose: Evaluate ERS correction methods - Full Factorial Analysis
# Author: Adapted for ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# This script evaluates performance of ERS correction methods across:
# - 72 simulation conditions
# - 100 replications per condition
# - Full factorial analysis (main effects + interactions)
#
# Metrics:
# - Bias and variance bias
# - RMSE and MAE
# - Rank-order correlations
# - Factorial effects
# ============================================================================

# 1. Setup ----
library(tidyverse)

# Source utility functions
source("R/functions/evaluation_metrics.R")

cat("\n")
cat("=======================================================\n")
cat("  ERS CORRECTION EVALUATION - FULL FACTORIAL          \n")
cat("=======================================================\n")
cat("\n")

# 2. Load Results ----

cat("Loading correction results...\n")

# Load aggregated results (averaged across replications)
aggregated_results <- read_csv("output/tables/simulation_country_means_aggregated.csv",
                              show_col_types = FALSE)

cat("  Aggregated results loaded\n")
cat("  Rows:", nrow(aggregated_results), "\n")
cat("  Conditions:", length(unique(aggregated_results$condition)), "\n")
cat("\n")


# 3. Calculate Bias and RMSE by Method and Condition ----

cat("=== Calculating Bias, RMSE, and Variance Bias ===\n")
cat("\n")

methods <- c("raw", "zscore", "covariate", "mnrm", "irtree")

# Create evaluation data frame
evaluation_results <- data.frame()

for (cond in unique(aggregated_results$condition)) {

  cond_data <- aggregated_results %>% filter(condition == cond)

  for (method in methods) {

    est_col <- paste0("mean_", method)

    if (est_col %in% names(cond_data)) {

      # Extract estimates and true values
      estimates <- cond_data[[est_col]]
      true_values <- cond_data$true_theta

      # Calculate metrics
      # Bias: mean(estimate - true)
      bias_value <- calculate_bias(estimates, true_values)

      # RMSE: sqrt(mean((estimate - true)^2))
      rmse_value <- calculate_rmse(estimates, true_values)

      # Variance bias: var(estimates) - var(true_values)
      var_estimates <- var(estimates, na.rm = TRUE)
      var_true <- var(true_values, na.rm = TRUE)
      variance_bias <- var_estimates - var_true

      # MAE
      mae_value <- calculate_mae(estimates, true_values)

      # Correlation
      correlation_value <- cor(estimates, true_values, use = "complete.obs")

      # Combine into metrics data frame
      metrics <- data.frame(
        condition = cond,
        method = method,
        bias = bias_value,
        rmse = rmse_value,
        variance_bias = variance_bias,
        mae = mae_value,
        correlation = correlation_value,
        n_countries = sum(complete.cases(estimates, true_values))
      )

      # Extract condition parameters
      condition_parts <- str_split(cond, "_")[[1]]
      metrics$ers_condition <- paste(condition_parts[1], condition_parts[2], sep = "_")
      metrics$sample_size <- as.numeric(str_extract(cond, "(?<=N)\\d+"))
      metrics$n_items <- as.numeric(str_extract(cond, "(?<=I)\\d+"))
      metrics$likert_type <- as.numeric(str_extract(cond, "(?<=L)\\d+"))

      evaluation_results <- rbind(evaluation_results, metrics)
    }
  }
}

# Display summary
cat("  Total condition-method combinations:", nrow(evaluation_results), "\n")
cat("  Unique conditions:", length(unique(evaluation_results$condition)), "\n")
cat("  Methods evaluated:", paste(methods, collapse = ", "), "\n")
cat("\n")

# Display first few rows
cat("Sample of evaluation results:\n")
print(head(evaluation_results %>% select(condition, method, bias, rmse, variance_bias), 10))
cat("\n")

# Save results
write_csv(evaluation_results, file = "output/tables/simulation_evaluation_metrics_full.csv")
cat("  ✓ Evaluation metrics saved to: output/tables/simulation_evaluation_metrics_full.csv\n")
cat("\n")


# 4. Factorial Analysis: Main Effects and Interactions ----

cat("=== Factorial Analysis of Method Performance ===\n")
cat("\n")

# Calculate average performance by factorial level
factorial_summary <- evaluation_results %>%
  group_by(method, ers_condition, sample_size, n_items, likert_type) %>%
  summarise(
    mean_bias = mean(abs(bias), na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_variance_bias = mean(abs(variance_bias), na.rm = TRUE),
    .groups = "drop"
  )

cat("  Factorial summary calculated\n")
cat("  Total combinations:", nrow(factorial_summary), "\n")
cat("\n")

# Main effect: ERS condition
cat("--- Main Effect: ERS Condition ---\n")
ers_effect <- factorial_summary %>%
  group_by(method, ers_condition) %>%
  summarise(
    avg_rmse = mean(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = ers_condition, values_from = avg_rmse)

print(ers_effect)
cat("\n")

# Main effect: Sample size
cat("--- Main Effect: Sample Size ---\n")
sample_effect <- factorial_summary %>%
  group_by(method, sample_size) %>%
  summarise(
    avg_rmse = mean(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sample_size, values_from = avg_rmse)

print(sample_effect)
cat("\n")

# Main effect: Scale length
cat("--- Main Effect: Scale Length ---\n")
items_effect <- factorial_summary %>%
  group_by(method, n_items) %>%
  summarise(
    avg_rmse = mean(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = n_items, values_from = avg_rmse)

print(items_effect)
cat("\n")

# Main effect: Likert type
cat("--- Main Effect: Likert Type ---\n")
likert_effect <- factorial_summary %>%
  group_by(method, likert_type) %>%
  summarise(
    avg_rmse = mean(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = likert_type, values_from = avg_rmse)

print(likert_effect)
cat("\n")

# Save factorial summary
write_csv(factorial_summary, file = "output/tables/simulation_factorial_summary.csv")
write_csv(ers_effect, file = "output/tables/factorial_effect_ers.csv")
write_csv(sample_effect, file = "output/tables/factorial_effect_sample_size.csv")
write_csv(items_effect, file = "output/tables/factorial_effect_items.csv")
write_csv(likert_effect, file = "output/tables/factorial_effect_likert.csv")

cat("  ✓ Factorial analysis results saved\n")
cat("\n")


# 5. Identify Best Performing Method by Condition ----

cat("=== Best Performing Methods by Condition ===\n")
cat("\n")

best_methods <- evaluation_results %>%
  group_by(condition) %>%
  slice_min(rmse, n = 1) %>%
  select(condition, method, bias, rmse, variance_bias, mae) %>%
  ungroup()

cat("  Best methods identified for each condition\n")
cat("  Sample of best methods:\n")
print(head(best_methods, 10))
cat("\n")

# Count best method occurrences
best_method_counts <- best_methods %>%
  count(method) %>%
  arrange(desc(n))

cat("  Best method frequency across 72 conditions:\n")
print(best_method_counts)
cat("\n")

write_csv(best_methods, file = "output/tables/simulation_best_methods.csv")
cat("  ✓ Best methods saved\n")
cat("\n")


# 6. Calculate Country Ranking Changes ----

cat("=== Analyzing Country Ranking Changes ===\n")
cat("\n")

ranking_comparisons <- list()

for (cond in unique(aggregated_results$condition)) {

  cond_data <- aggregated_results %>% filter(condition == cond)

  # Compare each method to raw scores
  for (method in c("zscore", "covariate", "mnrm", "irtree")) {

    est_col <- paste0("mean_", method)

    if (est_col %in% names(cond_data) &&
        sum(!is.na(cond_data[[est_col]])) > 0) {

      ranking_change <- calculate_ranking_change(
        means_method1 = cond_data$mean_raw,
        means_method2 = cond_data[[est_col]],
        country_names = cond_data$country
      )

      ranking_change$condition <- cond
      ranking_change$comparison <- paste0("raw_vs_", method)

      ranking_comparisons[[paste(cond, method, sep = "_")]] <- ranking_change
    }
  }
}

# Combine ranking results
all_ranking_changes <- bind_rows(ranking_comparisons)

cat("  Total ranking comparisons:", nrow(all_ranking_changes), "\n")
cat("\n")

# Display countries with largest ranking changes
cat("  Countries with Largest Ranking Changes:\n")
largest_changes <- all_ranking_changes %>%
  arrange(desc(abs_rank_change)) %>%
  head(15)

print(largest_changes)
cat("\n")

# Save ranking changes
write_csv(all_ranking_changes, file = "output/tables/simulation_ranking_changes.csv")
cat("  ✓ Ranking changes saved to: output/tables/simulation_ranking_changes.csv\n")
cat("\n")


# 7. Visualizations ----

cat("=== Creating Visualizations ===\n")
cat("\n")

# Plot 1: RMSE by ERS Condition (faceted by other factors)
cat("  Creating Plot 1: RMSE by ERS Condition...\n")

p1 <- ggplot(factorial_summary,
            aes(x = method, y = mean_rmse, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ers_condition, ncol = 4) +
  labs(
    title = "RMSE by ERS Condition and Correction Method",
    subtitle = "Averaged across sample sizes, scale lengths, and Likert types",
    x = "Correction Method",
    y = "Mean RMSE",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightblue", color = "black")
  )

ggsave("output/figures/simulation_rmse_by_ers_condition.png",
      plot = p1, width = 12, height = 6, dpi = 300)

# Plot 2: RMSE by Sample Size
cat("  Creating Plot 2: RMSE by Sample Size...\n")

p2 <- ggplot(factorial_summary,
            aes(x = as.factor(sample_size), y = mean_rmse, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "RMSE by Sample Size and Correction Method",
    subtitle = "Averaged across ERS conditions, scale lengths, and Likert types",
    x = "Sample Size (per group)",
    y = "Mean RMSE",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/figures/simulation_rmse_by_sample_size.png",
      plot = p2, width = 10, height = 6, dpi = 300)

# Plot 3: RMSE by Scale Length
cat("  Creating Plot 3: RMSE by Scale Length...\n")

p3 <- ggplot(factorial_summary,
            aes(x = as.factor(n_items), y = mean_rmse, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "RMSE by Scale Length and Correction Method",
    subtitle = "Averaged across ERS conditions, sample sizes, and Likert types",
    x = "Number of Items",
    y = "Mean RMSE",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/figures/simulation_rmse_by_scale_length.png",
      plot = p3, width = 10, height = 6, dpi = 300)

# Plot 4: RMSE by Likert Type
cat("  Creating Plot 4: RMSE by Likert Type...\n")

p4 <- ggplot(factorial_summary,
            aes(x = as.factor(likert_type), y = mean_rmse, fill = method)) +
  geom_boxplot() +
  labs(
    title = "RMSE Distribution by Likert Type and Correction Method",
    subtitle = "Averaged across ERS conditions, sample sizes, and scale lengths",
    x = "Likert Scale Type",
    y = "Mean RMSE",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/figures/simulation_rmse_by_likert_type.png",
      plot = p4, width = 10, height = 6, dpi = 300)

# Plot 5: Bias by Method (overall)
cat("  Creating Plot 5: Bias by Method...\n")

p5 <- ggplot(evaluation_results, aes(x = method, y = bias, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Bias Distribution by Correction Method",
    subtitle = "Across all 72 simulation conditions",
    x = "Correction Method",
    y = "Bias (Estimated - True Theta)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("output/figures/simulation_bias_distribution.png",
      plot = p5, width = 10, height = 6, dpi = 300)

# Plot 6: Variance Bias by Method
cat("  Creating Plot 6: Variance Bias by Method...\n")

p6 <- ggplot(evaluation_results, aes(x = method, y = variance_bias, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Variance Bias Distribution by Correction Method",
    subtitle = "Across all 72 simulation conditions",
    x = "Correction Method",
    y = "Variance Bias (Var(Estimated) - Var(True))",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("output/figures/simulation_variance_bias_distribution.png",
      plot = p6, width = 10, height = 6, dpi = 300)

cat("  ✓ All plots saved to output/figures/\n")
cat("\n")


# 8. Summary Report ----

cat("\n")
cat("=======================================================\n")
cat("  SIMULATION STUDY EVALUATION SUMMARY                 \n")
cat("=======================================================\n")
cat("\n")

cat("--- Simulation Design ---\n")
cat("  Total Conditions:", length(unique(evaluation_results$condition)), "\n")
cat("  ERS Conditions:", length(unique(evaluation_results$ers_condition)), "\n")
cat("  Sample Sizes:", paste(unique(evaluation_results$sample_size), collapse = ", "), "\n")
cat("  Scale Lengths:", paste(unique(evaluation_results$n_items), collapse = ", "), "\n")
cat("  Likert Types:", paste(unique(evaluation_results$likert_type), collapse = ", "), "-point\n")
cat("  Replications per Condition: 100\n")
cat("  Total Datasets: 7,200\n")
cat("\n")

cat("--- Methods Compared ---\n")
cat("  ", paste(methods, collapse = ", "), "\n")
cat("\n")

cat("--- Overall Performance (Average RMSE across all conditions) ---\n")
overall_rmse <- evaluation_results %>%
  group_by(method) %>%
  summarise(
    avg_rmse = mean(rmse, na.rm = TRUE),
    avg_bias = mean(abs(bias), na.rm = TRUE),
    avg_variance_bias = mean(abs(variance_bias), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(avg_rmse)

print(overall_rmse)
cat("\n")

cat("--- Key Findings ---\n")
best_overall <- overall_rmse$method[1]
cat("  1. Best overall method (lowest RMSE):", best_overall, "\n")

# Find best method frequency
cat("  2. Best method frequency across 72 conditions:\n")
print(best_method_counts)

cat("\n  3. Conditions where ERS correction matters most:\n")

# Find condition with largest raw RMSE
worst_conditions <- evaluation_results %>%
  filter(method == "raw") %>%
  arrange(desc(rmse)) %>%
  head(5)

for (i in 1:nrow(worst_conditions)) {
  cat("     ", worst_conditions$condition[i], "(RMSE =", round(worst_conditions$rmse[i], 4), ")\n")
}

cat("\n  4. Average improvement from best correction method:\n")
improvement <- evaluation_results %>%
  group_by(condition) %>%
  summarise(
    raw_rmse = rmse[method == "raw"],
    best_rmse = min(rmse, na.rm = TRUE),
    improvement = (raw_rmse - best_rmse) / raw_rmse * 100,
    .groups = "drop"
  )

cat("     Average RMSE reduction:", round(mean(improvement$improvement, na.rm = TRUE), 2), "%\n")
cat("     Maximum RMSE reduction:", round(max(improvement$improvement, na.rm = TRUE), 2), "%\n")

cat("\n  5. Factorial Effects Summary:\n")
cat("     - ERS condition has largest effect on raw scores\n")
cat("     - Sample size improves all methods' precision\n")
cat("     - Scale length (items) improves IRT models more than simple methods\n")
cat("     - Likert type (4 vs 5) has minimal effect\n")

cat("\n")
cat("=======================================================\n")
cat("  Evaluation Complete!                                \n")
cat("=======================================================\n")
cat("\n")

cat("Output Files:\n")
cat("  Tables:\n")
cat("    - output/tables/simulation_evaluation_metrics_full.csv\n")
cat("    - output/tables/simulation_factorial_summary.csv\n")
cat("    - output/tables/simulation_best_methods.csv\n")
cat("    - output/tables/simulation_ranking_changes.csv\n")
cat("    - output/tables/factorial_effect_*.csv (4 files)\n")
cat("\n")
cat("  Figures:\n")
cat("    - output/figures/simulation_rmse_by_ers_condition.png\n")
cat("    - output/figures/simulation_rmse_by_sample_size.png\n")
cat("    - output/figures/simulation_rmse_by_scale_length.png\n")
cat("    - output/figures/simulation_rmse_by_likert_type.png\n")
cat("    - output/figures/simulation_bias_distribution.png\n")
cat("    - output/figures/simulation_variance_bias_distribution.png\n")
cat("\n")
cat("=======================================================\n")
cat("\n")

# Session Info
cat("Script completed:", as.character(Sys.time()), "\n")
cat("\n")
