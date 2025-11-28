# ============================================================================
# Script: 03_evaluation.R
# Purpose: Evaluate ERS correction methods (Bias, RMSE, ranking changes)
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)

# Source utility functions
source("R/functions/evaluation_metrics.R")

# 2. Load Results ----

cat("Loading correction results...\n")
all_results <- readRDS("data/simulated/correction_results.rds")


# 3. Calculate Bias and RMSE by Method ----

cat("\n=== Calculating Bias and RMSE for Each Method ===\n")

# Combine all conditions
combined_data <- bind_rows(all_results, .id = "condition")

# Create evaluation data frame
evaluation_results <- data.frame()

methods <- c("raw", "zscore", "covariate", "mnrm", "irtree")

for (cond in unique(combined_data$condition)) {

  cond_data <- combined_data %>% filter(condition == cond)

  for (method in methods) {

    est_col <- paste0("mean_", method)

    if (est_col %in% names(cond_data)) {

      # Extract estimates and true values
      estimates <- cond_data[[est_col]]
      true_values <- cond_data$true_theta

      # Calculate metrics
      metrics <- data.frame(
        condition = cond,
        method = method,
        bias = calculate_bias(estimates, true_values),
        rmse = calculate_rmse(estimates, true_values),
        mae = calculate_mae(estimates, true_values),
        correlation = cor(estimates, true_values, use = "complete.obs"),
        n_countries = sum(complete.cases(estimates, true_values))
      )

      evaluation_results <- rbind(evaluation_results, metrics)
    }
  }
}


# Display results
cat("\nBias and RMSE by Condition and Method:\n")
print(evaluation_results %>% arrange(condition, rmse))

# Save results
write_csv(evaluation_results,
         file = "output/tables/simulation_evaluation_metrics.csv")

cat("\nEvaluation metrics saved to output/tables/simulation_evaluation_metrics.csv\n")


# 4. Identify Best Performing Method by Condition ----

cat("\n=== Best Performing Methods by Condition ===\n")

best_methods <- evaluation_results %>%
  group_by(condition) %>%
  slice_min(rmse, n = 1) %>%
  select(condition, method, bias, rmse, mae) %>%
  ungroup()

print(best_methods)


# 5. Calculate Country Ranking Changes ----

cat("\n=== Analyzing Country Ranking Changes ===\n")

ranking_comparisons <- list()

for (cond in unique(combined_data$condition)) {

  cond_data <- combined_data %>% filter(condition == cond)

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

# Display countries with largest ranking changes
cat("\nCountries with Largest Ranking Changes:\n")
largest_changes <- all_ranking_changes %>%
  arrange(desc(abs_rank_change)) %>%
  head(20)

print(largest_changes)

# Save ranking changes
write_csv(all_ranking_changes,
         file = "output/tables/simulation_ranking_changes.csv")

cat("\nRanking changes saved to output/tables/simulation_ranking_changes.csv\n")


# 6. Visualize Bias and RMSE ----

cat("\n=== Creating visualization of evaluation metrics ===\n")

# Create plot
p1 <- ggplot(evaluation_results, aes(x = method, y = rmse, fill = method)) +
  geom_col() +
  facet_wrap(~condition, scales = "free_y") +
  labs(
    title = "RMSE by ERS Correction Method and Simulation Condition",
    x = "Correction Method",
    y = "Root Mean Square Error (RMSE)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("output/figures/simulation_rmse_comparison.png",
      plot = p1, width = 10, height = 6, dpi = 300)

cat("Plot saved to output/figures/simulation_rmse_comparison.png\n")


# Create bias plot
p2 <- ggplot(evaluation_results, aes(x = method, y = bias, fill = method)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~condition, scales = "free_y") +
  labs(
    title = "Bias by ERS Correction Method and Simulation Condition",
    x = "Correction Method",
    y = "Bias (Estimated - True)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("output/figures/simulation_bias_comparison.png",
      plot = p2, width = 10, height = 6, dpi = 300)

cat("Plot saved to output/figures/simulation_bias_comparison.png\n")


# 7. Summary Report ----

cat("\n\n")
cat("=======================================================\n")
cat("        SIMULATION STUDY EVALUATION SUMMARY           \n")
cat("=======================================================\n\n")

cat("Total Conditions Evaluated:", length(unique(combined_data$condition)), "\n")
cat("Methods Compared:", paste(methods, collapse = ", "), "\n")
cat("Countries per Condition:", unique(combined_data %>% count(condition))$n, "\n\n")

cat("--- Overall Performance (Average RMSE across conditions) ---\n")
overall_rmse <- evaluation_results %>%
  group_by(method) %>%
  summarise(
    avg_rmse = mean(rmse, na.rm = TRUE),
    avg_bias = mean(abs(bias), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(avg_rmse)

print(overall_rmse)

cat("\n--- Key Findings ---\n")
best_overall <- overall_rmse$method[1]
cat("1. Best overall method:", best_overall, "\n")
cat("2. Conditions where ERS correction matters most:\n")

# Find condition with largest raw RMSE
worst_condition <- evaluation_results %>%
  filter(method == "raw") %>%
  arrange(desc(rmse)) %>%
  head(1)

cat("   ", worst_condition$condition, "(RMSE =", round(worst_condition$rmse, 4), ")\n")

cat("\n3. Average improvement from best correction method:\n")
improvement <- evaluation_results %>%
  group_by(condition) %>%
  summarise(
    raw_rmse = rmse[method == "raw"],
    best_rmse = min(rmse),
    improvement = (raw_rmse - best_rmse) / raw_rmse * 100,
    .groups = "drop"
  )

cat("   Average RMSE reduction:", round(mean(improvement$improvement, na.rm = TRUE), 2), "%\n")

cat("\n=======================================================\n")
cat("Evaluation complete. Results saved to output/tables/\n")
cat("Figures saved to output/figures/\n")
cat("=======================================================\n\n")

# 8. Session Info ----
cat("Script completed:", as.character(Sys.time()), "\n")
