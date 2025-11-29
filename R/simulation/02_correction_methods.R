# ============================================================================
# Script: 02_correction_methods.R
# Purpose: Apply ERS correction methods to simulated data (Full Factorial)
# Author: Adapted for ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# This script processes all 72 simulation conditions × 100 replications
# with varying:
# - ERS conditions (4)
# - Sample sizes (100, 500, 1000)
# - Scale lengths (10, 20, 30 items)
# - Likert types (4-point, 5-point)
#
# Total: 7,200 datasets to process
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(mirt)

# Source utility functions
source("R/functions/ers_indices.R")
source("R/functions/correction_methods.R")
source("R/functions/irtree_functions.R")

set.seed(12345)

# Initialize IRTree custom nodes
irtree_nodes <- initialize_irtree_nodes()

cat("\n")
cat("=======================================================\n")
cat("  ERS CORRECTION METHODS - FULL FACTORIAL             \n")
cat("=======================================================\n")
cat("\n")

# 2. Load Simulated Data ----

cat("Loading simulated datasets...\n")
all_datasets <- readRDS("data/simulated/all_conditions_full_factorial.rds")

cat("  Total conditions:", length(all_datasets), "\n")
cat("  Replications per condition:", length(all_datasets[[1]]), "\n")
cat("\n")


# 3. Helper Functions ----

#' Extract Condition Parameters from Condition Name
#'
#' @param condition_name String like "low_ers_N500_I10_L4"
#' @return List with ers_condition, sample_size, n_items, likert_type
parse_condition_name <- function(condition_name) {

  # Extract components using regex
  ers_condition <- str_extract(condition_name, "^[^_]+_[^_]+")
  sample_size <- as.numeric(str_extract(condition_name, "(?<=_N)\\d+"))
  n_items <- as.numeric(str_extract(condition_name, "(?<=_I)\\d+"))
  likert_type <- as.numeric(str_extract(condition_name, "(?<=_L)\\d+"))

  return(list(
    ers_condition = ers_condition,
    sample_size = sample_size,
    n_items = n_items,
    likert_type = likert_type
  ))
}


#' Apply All Correction Methods to a Dataset
#'
#' @param data Data frame with country, true_theta, true_ers, and Item columns
#' @param n_items Number of items in this dataset
#' @param n_categories Number of response categories (4 or 5)
#' @return Data frame with country means for each method
apply_all_corrections <- function(data, n_items, n_categories) {

  # Determine item columns based on n_items
  item_cols <- paste0("Item", 1:n_items)

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
  model_syntax <- paste0(
    "Theta = 1-", n_items, "\n",
    "ERS = 1-", n_items, "\n",
    "FREE = (GROUP, COV_21)"
  )

  # Generate s matrix for MNRM (depends on number of categories)
  s_content <- matrix(0:(n_categories - 1), nrow = 1)
  s_ers <- matrix(rep(0, n_categories), nrow = 1)
  s_ers[1, 1] <- 1  # First category (strongly disagree) = extreme
  s_ers[1, n_categories] <- 1  # Last category (strongly agree) = extreme
  s_matrix <- rbind(s_content, s_ers)

  s_mirt_list <- vector("list", length = n_items)
  for (j in 1:n_items) {
    s_mirt_list[[j]] <- t(s_matrix)
  }

  # Fit MNRM model (may take time)
  cat("    Fitting MNRM model (", n_items, " items, ", n_categories, " categories)...\n", sep = "")

  # Initialize model3_means in case of errors
  model3_means <- data.frame(
    country = unique(group_var),
    mean_mnrm = NA
  )

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


  # === Model 4: IRTree ===
  # Note: IRTree is the most complex model, computationally intensive
  # Uses irtrees package approach with pseudo-item decomposition

  cat("    Fitting IRTree model...\n")

  # Initialize model4_means in case of errors
  model4_means <- data.frame(
    country = unique(group_var),
    mean_irtree = NA
  )

  tryCatch({
    # For IRTree, we use a simplified multidimensional approach
    # Full IRTree requires pseudo-item decomposition (see irtrees package)
    # Here we approximate with structured GPCM

    # IRTree-like model with constrained parameters
    irtree_syntax <- paste0(
      "Theta = 1-", n_items, "\n",
      "ERS = 1-", n_items, "\n",
      "FREE = (GROUP, COV_21)"
    )

    # Use graded response model as approximation
    # (Full IRTree would require custom node implementation)
    irtree_fit <- multipleGroup(
      data = mirt_data,
      model = irtree_syntax,
      group = group_var,
      itemtype = "graded",  # Approximation
      method = "EM",
      invariance = c("free_mean", "free_var"),
      technical = list(NCYCLES = 1000),
      verbose = FALSE
    )

    # Extract means
    coef_list_irtree <- coef(irtree_fit, simplify = TRUE)
    group_names <- unique(group_var)

    model4_estimates <- sapply(group_names, function(g_name) {
      group_idx <- which(names(coef_list_irtree) == g_name)
      if (length(group_idx) > 0 && "means" %in% names(coef_list_irtree[[group_idx]])) {
        return(coef_list_irtree[[group_idx]]$means[1])
      } else {
        return(NA)
      }
    })

    model4_means <- data.frame(
      country = group_names,
      mean_irtree = model4_estimates
    )

  }, error = function(e) {
    cat("    IRTree fitting failed:", e$message, "\n")
    cat("    Using NA for IRTree estimates\n")
    cat("    Note: Full IRTree implementation requires irtrees package\n")
    model4_means <- data.frame(
      country = unique(group_var),
      mean_irtree = NA
    )
  })


  # === Combine All Results ===
  results <- model0_means %>%
    left_join(model1_means, by = "country") %>%
    left_join(model2_means, by = "country") %>%
    left_join(model3_means, by = "country") %>%
    left_join(model4_means, by = "country")

  return(results)
}


# 4. Apply Corrections to All Conditions and Replications ----

cat("Applying ERS correction methods to all conditions and replications...\n")
cat("This will process 7,200 datasets - please be patient!\n")
cat("\n")

all_results <- list()
condition_counter <- 0
start_time_total <- Sys.time()

for (condition_name in names(all_datasets)) {

  condition_counter <- condition_counter + 1

  cat("\n")
  cat("--- Condition", condition_counter, "/", length(all_datasets), ":", condition_name, "---\n")

  # Parse condition parameters
  condition_params <- parse_condition_name(condition_name)
  n_items <- condition_params$n_items
  n_categories <- condition_params$likert_type

  cat("  ERS condition:", condition_params$ers_condition, "\n")
  cat("  Sample size:", condition_params$sample_size, "per group\n")
  cat("  Items:", n_items, "\n")
  cat("  Likert scale:", n_categories, "-point\n")
  cat("  Processing", length(all_datasets[[condition_name]]), "replications...\n")

  start_time_condition <- Sys.time()

  # Get all replications for this condition
  replications <- all_datasets[[condition_name]]
  n_replications <- length(replications)

  # Storage for this condition's results
  condition_results <- list()

  # Process each replication
  for (rep in 1:n_replications) {

    # Progress indicator (every 25 reps)
    if (rep %% 25 == 0) {
      cat("    Rep", rep, "/", n_replications, "\n")
    }

    # Get data for this replication
    sim_data <- replications[[rep]]$data

    # Apply all correction methods
    results <- apply_all_corrections(sim_data, n_items, n_categories)

    # Add metadata
    results$condition <- condition_name
    results$replication <- rep
    results$ers_condition <- condition_params$ers_condition
    results$sample_size <- condition_params$sample_size
    results$n_items <- n_items
    results$likert_type <- n_categories

    # Store
    condition_results[[rep]] <- results
  }

  # Combine replications for this condition
  all_results[[condition_name]] <- bind_rows(condition_results)

  # Report timing
  time_elapsed <- difftime(Sys.time(), start_time_condition, units = "secs")
  cat("  ✓ Completed in", round(time_elapsed, 1), "seconds\n")
}

# Combine all results
cat("\nCombining all results...\n")
combined_results <- bind_rows(all_results)

cat("  Total rows:", nrow(combined_results), "\n")
cat("  Total conditions:", length(unique(combined_results$condition)), "\n")
cat("  Total replications:", max(combined_results$replication), "\n")


# 5. Save Results ----

cat("\n")
cat("=======================================================\n")
cat("  Saving Results                                      \n")
cat("=======================================================\n")
cat("\n")

# Save full results (by condition)
saveRDS(all_results, file = "data/simulated/correction_results_full_factorial.rds")
cat("  ✓ Full results saved to: data/simulated/correction_results_full_factorial.rds\n")

# Save combined results (all replications)
write_csv(combined_results, file = "output/tables/simulation_country_means_full.csv")
cat("  ✓ Combined results saved to: output/tables/simulation_country_means_full.csv\n")


# 6. Aggregate Across Replications ----

cat("\n")
cat("=======================================================\n")
cat("  Aggregating Across Replications                    \n")
cat("=======================================================\n")
cat("\n")

# Calculate mean estimates across replications for each condition
aggregated_results <- combined_results %>%
  group_by(condition, ers_condition, sample_size, n_items, likert_type, country) %>%
  summarise(
    # Average estimates across 100 replications
    mean_raw = mean(mean_raw, na.rm = TRUE),
    mean_zscore = mean(mean_zscore, na.rm = TRUE),
    mean_covariate = mean(mean_covariate, na.rm = TRUE),
    mean_mnrm = mean(mean_mnrm, na.rm = TRUE),
    mean_irtree = mean(mean_irtree, na.rm = TRUE),
    # True values (should be constant across reps)
    true_theta = mean(true_theta, na.rm = TRUE),
    true_ers = mean(true_ers, na.rm = TRUE),
    # Number of valid replications
    n_reps = n(),
    .groups = "drop"
  )

cat("  Aggregated results by condition:\n")
cat("    Total condition-country combinations:", nrow(aggregated_results), "\n")

# Save aggregated results
write_csv(aggregated_results, file = "output/tables/simulation_country_means_aggregated.csv")
cat("  ✓ Aggregated results saved to: output/tables/simulation_country_means_aggregated.csv\n")


# 7. Summary Comparison ----

cat("\n")
cat("=======================================================\n")
cat("  Summary: Method Performance                        \n")
cat("=======================================================\n")
cat("\n")

# Reshape for easier comparison (using aggregated data)
comparison_long <- aggregated_results %>%
  pivot_longer(
    cols = c(mean_raw, mean_zscore, mean_covariate, mean_mnrm, mean_irtree),
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
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    .groups = "drop"
  )

cat("  Summary statistics calculated for all methods\n")
cat("  Total condition-method combinations:", nrow(comparison_summary), "\n")

write_csv(comparison_summary, file = "output/tables/simulation_method_comparison.csv")
cat("  ✓ Method comparison saved to: output/tables/simulation_method_comparison.csv\n")


# 8. Final Summary ----

total_time <- difftime(Sys.time(), start_time_total, units = "mins")

cat("\n")
cat("=======================================================\n")
cat("  ERS Correction Methods Completed!                   \n")
cat("=======================================================\n")
cat("\n")
cat("Total conditions processed:", length(all_datasets), "\n")
cat("Total replications per condition:", 100, "\n")
cat("Total datasets processed:", nrow(combined_results) / 2, "\n")  # Divide by 2 countries
cat("Total time:", round(total_time, 2), "minutes\n")
cat("Average time per condition:", round(total_time / length(all_datasets), 2), "minutes\n")
cat("\n")

cat("Next step:\n")
cat("  Run: source('R/simulation/03_evaluation.R')\n")
cat("\n")
cat("=======================================================\n")
