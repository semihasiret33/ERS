# ============================================================================
# Script: country_comparisons.R
# Purpose: Create visualizations comparing country means across ERS methods
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(ggplot2)
library(scales)

# 2. Load Results ----

cat("Loading PISA analysis results...\n")
all_models_results <- readRDS("data/processed/pisa_all_models_results.rds")


# 3. Country Means Comparison Plot ----

cat("Creating country means comparison plot...\n")

# Reshape data to long format for plotting
plot_data <- all_models_results %>%
  select(CNT, mean_model0, mean_model1, mean_model2, mean_model3, mean_model4) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "model",
    values_to = "country_mean"
  ) %>%
  mutate(
    model = factor(model,
                  levels = c("mean_model0", "mean_model1",
                            "mean_model2", "mean_model3", "mean_model4"),
                  labels = c("Model 0: Raw", "Model 1: Z-score",
                            "Model 2: Covariate", "Model 3: ML-MNRM", "Model 4: IRTree"))
  )

# Create plot
p1 <- ggplot(plot_data, aes(x = CNT, y = country_mean, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Country Means Across ERS Correction Methods",
    subtitle = "Comparison of raw scores and ERS-corrected estimates",
    x = "Country",
    y = "Mean Score",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("output/figures/country_means_comparison.png",
      plot = p1, width = 10, height = 6, dpi = 300)

cat("Saved: output/figures/country_means_comparison.png\n")


# 4. Ranking Changes Plot ----

cat("Creating ranking changes plot...\n")

# Reshape ranking data
rank_data <- all_models_results %>%
  select(CNT, rank_model0, rank_model1, rank_model2, rank_model3, rank_model4) %>%
  pivot_longer(
    cols = starts_with("rank_"),
    names_to = "model",
    values_to = "rank"
  ) %>%
  mutate(
    model = factor(model,
                  levels = c("rank_model0", "rank_model1",
                            "rank_model2", "rank_model3", "rank_model4"),
                  labels = c("Model 0", "Model 1", "Model 2", "Model 3", "Model 4"))
  )

# Create line plot showing ranking changes
p2 <- ggplot(rank_data, aes(x = model, y = rank, group = CNT, color = CNT)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_y_reverse(breaks = 1:nrow(all_models_results)) +
  labs(
    title = "Country Ranking Changes Across ERS Correction Methods",
    subtitle = "Lines show how each country's rank changes with different corrections",
    x = "Model",
    y = "Rank (1 = Highest)",
    color = "Country"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("output/figures/ranking_changes.png",
      plot = p2, width = 10, height = 6, dpi = 300)

cat("Saved: output/figures/ranking_changes.png\n")


# 5. ERS vs Country Means Scatterplot ----

cat("Creating ERS vs country means scatterplot...\n")

p3 <- ggplot(all_models_results,
            aes(x = mean_ers_greenleaf, y = mean_model0, label = CNT)) +
  geom_point(size = 4, alpha = 0.7, color = "steelblue") +
  geom_text(vjust = -1, size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
  labs(
    title = "Relationship Between ERS Level and Country Mean Scores",
    subtitle = "Higher ERS may inflate or deflate observed country means",
    x = "ERS Level (Greenleaf Index)",
    y = "Country Mean (Raw Scores)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("output/figures/ers_vs_country_means.png",
      plot = p3, width = 8, height = 6, dpi = 300)

cat("Saved: output/figures/ers_vs_country_means.png\n")


# 6. Heatmap of Rank Changes ----

cat("Creating heatmap of rank changes...\n")

# Prepare data for heatmap
heatmap_data <- all_models_results %>%
  select(CNT, rank_change_m1, rank_change_m2, rank_change_m3, rank_change_m4) %>%
  pivot_longer(
    cols = starts_with("rank_change_"),
    names_to = "comparison",
    values_to = "rank_change"
  ) %>%
  mutate(
    comparison = factor(comparison,
                       levels = c("rank_change_m1", "rank_change_m2",
                                 "rank_change_m3", "rank_change_m4"),
                       labels = c("Model 0 vs 1", "Model 0 vs 2",
                                 "Model 0 vs 3", "Model 0 vs 4"))
  )

p4 <- ggplot(heatmap_data,
            aes(x = comparison, y = CNT, fill = rank_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(rank_change, 1)), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Rank Change"
  ) +
  labs(
    title = "Country Rank Changes Relative to Raw Scores",
    subtitle = "Positive values = lower rank (worse) after correction, Negative = higher rank (better)",
    x = "Model Comparison",
    y = "Country"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("output/figures/rank_changes_heatmap.png",
      plot = p4, width = 8, height = 6, dpi = 300)

cat("Saved: output/figures/rank_changes_heatmap.png\n")


# 7. Method Agreement Plot ----

cat("Creating method agreement plot...\n")

# Calculate pairwise correlations (adding Model 4 comparisons)
cors <- data.frame(
  model1 = c("Model 0", "Model 0", "Model 0", "Model 0", "Model 1", "Model 1", "Model 1", "Model 2", "Model 2", "Model 3"),
  model2 = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 2", "Model 3", "Model 4", "Model 3", "Model 4", "Model 4"),
  correlation = c(
    cor(all_models_results$mean_model0, all_models_results$mean_model1, use = "complete.obs"),
    cor(all_models_results$mean_model0, all_models_results$mean_model2, use = "complete.obs"),
    cor(all_models_results$mean_model0, all_models_results$mean_model3, use = "complete.obs"),
    cor(all_models_results$mean_model0, all_models_results$mean_model4, use = "complete.obs"),
    cor(all_models_results$mean_model1, all_models_results$mean_model2, use = "complete.obs"),
    cor(all_models_results$mean_model1, all_models_results$mean_model3, use = "complete.obs"),
    cor(all_models_results$mean_model1, all_models_results$mean_model4, use = "complete.obs"),
    cor(all_models_results$mean_model2, all_models_results$mean_model3, use = "complete.obs"),
    cor(all_models_results$mean_model2, all_models_results$mean_model4, use = "complete.obs"),
    cor(all_models_results$mean_model3, all_models_results$mean_model4, use = "complete.obs")
  )
)

p5 <- ggplot(cors, aes(x = paste(model1, "vs", model2),
                      y = correlation)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(correlation, 3)), vjust = -0.5, size = 4) +
  ylim(0, 1.05) +
  labs(
    title = "Agreement Between ERS Correction Methods",
    subtitle = "Spearman correlation of country means across methods",
    x = "Method Comparison",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("output/figures/method_agreement.png",
      plot = p5, width = 10, height = 6, dpi = 300)

cat("Saved: output/figures/method_agreement.png\n")


# 8. Summary ----

cat("\n=======================================================\n")
cat("All visualizations created successfully!\n")
cat("=======================================================\n\n")
cat("Saved figures:\n")
cat("  1. country_means_comparison.png (5 models)\n")
cat("  2. ranking_changes.png (5 models)\n")
cat("  3. ers_vs_country_means.png\n")
cat("  4. rank_changes_heatmap.png (4 comparisons)\n")
cat("  5. method_agreement.png (10 pairwise correlations)\n")
cat("\nAll figures saved to: output/figures/\n")
cat("=======================================================\n\n")

# 9. Session Info ----
cat("Script completed:", as.character(Sys.time()), "\n")
