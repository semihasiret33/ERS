# ============================================================================
# Script: 02_correction_methods.R
# Purpose: Apply ERS correction methods to simulated data
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(mirt)

# Source utility functions
source("R/functions/ers_indices.R")
source("R/functions/correction_methods.R")

set.seed(12345)

# 2. Load Simulated Data ----

cat("Loading simulated datasets...\n")
simulation_datasets <- readRDS("data/simulated/all_conditions.rds")

# Items to use for analysis
item_cols <- paste0("Item", 1:10)


# 3. Functions ----

#' Apply All Correction Methods to a Dataset
#'
#' @param data Data frame with country, true_theta, true_ers, and Item columns
#' @param item_cols Character vector of item column names
#' @return Data frame with country means for each method
apply_all_corrections <- function(data, item_cols) {

  # === Model 0: Raw Scores (Baseline) ===
  model0_means <- data %>%
    group_by(country) %>%
    summarise(
      mean_raw = mean(rowMeans(across(all_of(item_cols)), na.rm = TRUE)),
      true_theta = mean(true_theta),
      true_ers = mean(true_ers),
      .groups = "drop"
    )


  # === Model 1: Z-Score Standardization ===
  # Person-mean centering
  data_zscore <- apply_zscore_correction(
    data = data,
    items = item_cols,
    suffix = "_z"
  )

  # Calculate means of z-scored items
  zscore_cols <- paste0(item_cols, "_z")
  model1_means <- data_zscore %>%
    group_by(country) %>%
    summarise(
      mean_zscore = mean(rowMeans(across(all_of(zscore_cols)), na.rm = TRUE)),
      .groups = "drop"
    )


  # === Model 2: Covariate Control ===
  # Add ERS index
  data_ers <- add_ers_index(data, items = item_cols, index_name = "ers_index")

  # Calculate raw scale score
  data_ers$scale_score_raw <- rowMeans(data_ers[, item_cols], na.rm = TRUE)

  # Apply covariate correction
  data_covariate <- apply_covariate_correction(
    data = data_ers,
    items = item_cols,
    ers_index = "ers_index",
    method = "mean",
    return_residuals = TRUE
  )

  model2_means <- data_covariate %>%
    group_by(country) %>%
    summarise(
      mean_covariate = mean(scale_score_ers_corrected, na.rm = TRUE),
      .groups = "drop"
    )


  # === Model 3: MNRM (ML-MNRM) ===
  # Note: This is computationally intensive, so we'll use a simplified version
  # For full implementation, see analysis/06_model3_mlmnrm.R

  # Prepare data for mirt
  mirt_data <- data[, item_cols]
  group_var <- data$country

  # Define MNRM model
  # Two dimensions: content (Theta) and ERS
  n_items <- length(item_cols)
  model_syntax <- paste0(
    "Theta = 1-", n_items, "\n",
    "ERS = 1-", n_items, "\n",
    "FREE = (GROUP, COV_21)"
  )

  # Generate s matrix for MNRM
  s_content <- matrix(0:3, nrow = 1)  # 4 categories
  s_ers <- matrix(c(1, 0, 0, 1), nrow = 1)
  s_matrix <- rbind(s_content, s_ers)

  s_mirt_list <- vector("list", length = n_items)
  for (j in 1:n_items) {
    s_mirt_list[[j]] <- t(s_matrix)
  }

  # Fit MNRM model (may take time)
  cat("    Fitting MNRM model...\n")

  tryCatch({
    mnrm_fit <- multipleGroup(
      data = mirt_data,
      model = model_syntax,
      group = group_var,
      itemtype = "gpcm",
      gpcm_mats = s_mirt_list,
      method = "EM",
      invariance = c("free_mean", "free_var", "slopes", "intercepts"),
      technical = list(NCYCLES = 1000),
      verbose = FALSE
    )

    # Extract latent means by group
    coef_list <- coef(mnrm_fit, simplify = TRUE)
    group_names <- unique(group_var)

    model3_estimates <- sapply(group_names, function(g_name) {
      group_idx <- which(names(coef_list) == g_name)
      if (length(group_idx) > 0) {
        # Extract theta mean (first dimension)
        return(coef_list[[group_idx]]$means[1])
      } else {
        return(NA)
      }
    })

    model3_means <- data.frame(
      country = group_names,
      mean_mnrm = model3_estimates
    )

  }, error = function(e) {
    cat("    MNRM fitting failed:", e$message, "\n")
    cat("    Using NA for MNRM estimates\n")
    model3_means <- data.frame(
      country = unique(group_var),
      mean_mnrm = NA
    )
  })


  # === Combine All Results ===
  results <- model0_means %>%
    left_join(model1_means, by = "country") %>%
    left_join(model2_means, by = "country") %>%
    left_join(model3_means, by = "country")

  return(results)
}


# 4. Apply Corrections to All Conditions ----

cat("\nApplying ERS correction methods to all simulation conditions...\n")

all_results <- list()

for (condition_name in names(simulation_datasets)) {

  cat("\n--- Processing condition:", condition_name, "---\n")

  sim_data <- simulation_datasets[[condition_name]]$data

  # Apply all correction methods
  results <- apply_all_corrections(sim_data, item_cols)

  # Add condition name
  results$condition <- condition_name

  # Store
  all_results[[condition_name]] <- results

  # Display results
  cat("\nCountry means by method:\n")
  print(results)
}

# Combine all results
combined_results <- bind_rows(all_results)


# 5. Save Results ----

saveRDS(all_results, file = "data/simulated/correction_results.rds")
write_csv(combined_results, file = "output/tables/simulation_country_means.csv")

cat("\nResults saved to:\n")
cat("  - data/simulated/correction_results.rds\n")
cat("  - output/tables/simulation_country_means.csv\n")


# 6. Summary Comparison ----

cat("\n=== Summary: Country Mean Estimates Across Methods ===\n")

# Reshape for easier comparison
comparison_long <- combined_results %>%
  pivot_longer(
    cols = c(mean_raw, mean_zscore, mean_covariate, mean_mnrm),
    names_to = "method",
    values_to = "estimated_mean"
  ) %>%
  mutate(method = str_remove(method, "mean_"))

# Compare to true theta
comparison_summary <- comparison_long %>%
  mutate(
    error = estimated_mean - true_theta,
    abs_error = abs(error)
  ) %>%
  group_by(condition, method) %>%
  summarise(
    mean_bias = mean(error, na.rm = TRUE),
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

print(comparison_summary)

write_csv(comparison_summary,
         file = "output/tables/simulation_method_comparison.csv")


cat("\n=== Script completed ===\n")
