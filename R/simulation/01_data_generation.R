# ============================================================================
# Script: 01_data_generation.R
# Purpose: Generate simulation data for all conditions
# Author: ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# This script generates simulated data for all 12 simulation conditions:
# - 4 ERS conditions (low, medium, high_diff, equal)
# - 3 Sample sizes (N=100, N=500, N=1000)
# - 500 replications per condition
#
# Total: 12 conditions × 500 replications = 6,000 datasets
#
# Data generation uses MNRM-Falk approach (Falk & Cai, 2016) with:
# - Alpha = 1.5 for both dimensions
# - Varying item difficulties (m values)
# - Two threshold sets (average vs difficult)
# ============================================================================

# 1. Setup ----
library(MASS)
library(tidyverse)

# Source data generation functions
source("R/simulation/00_mnrm_falk_generator.R")

set.seed(12345)

cat("\n")
cat("=======================================================\n")
cat("  ERS SIMULATION DATA GENERATION                      \n")
cat("=======================================================\n")
cat("\n")


# 2. Define Simulation Parameters ----

cat("Simulation parameters:\n")
cat("  - Number of items: 10\n")
cat("  - Categories: 4 (Strongly Disagree to Strongly Agree)\n")
cat("  - Alpha (discrimination): 1.5 for both theta and ERS\n")
cat("  - Thresholds: [-1, 0, 1] (average difficulty)\n")
cat("  - Replications per condition: 500\n")
cat("  - Sample sizes: 100, 500, 1000\n")
cat("  - ERS conditions: 4\n")
cat("  - Total conditions: 12\n")
cat("  - Total datasets: 6,000\n")
cat("\n")

# Simulation parameters
N_ITEMS <- 10
CATEGORIES <- 4
ALPHA <- c(1.5, 1.5)
THRESHOLDS_AVERAGE <- c(-1, 0, 1)
N_REPLICATIONS <- 500
SAMPLE_SIZES <- c(100, 500, 1000)


# 3. Define ERS Conditions ----

cat("ERS Conditions:\n")
cat("\n")

# Condition 1: Low ERS (both groups low, no difference)
cat("  1. low_ers:\n")
cat("     - Group 1: theta~N(0,1), ERS~N(0,1)\n")
cat("     - Group 2: theta~N(0,1), ERS~N(0,1)\n")
cat("     - Purpose: Baseline with minimal ERS, no group difference\n")
cat("\n")

low_ers <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 0),
  ers_sd = c(1, 1)
)

# Condition 2: Medium ERS (group difference of 0.5 SD)
cat("  2. medium_ers:\n")
cat("     - Group 1: theta~N(0,1), ERS~N(0,1)\n")
cat("     - Group 2: theta~N(0,1), ERS~N(0.5,1)\n")
cat("     - Purpose: Moderate ERS difference between groups\n")
cat("\n")

medium_ers <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 0.5),
  ers_sd = c(1, 1)
)

# Condition 3: High ERS difference (1 SD difference)
cat("  3. high_ers_diff:\n")
cat("     - Group 1: theta~N(0,1), ERS~N(0,1)\n")
cat("     - Group 2: theta~N(0,1), ERS~N(1.0,1)\n")
cat("     - Purpose: Large ERS difference, theta equal\n")
cat("\n")

high_ers_diff <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 1.0),
  ers_sd = c(1, 1)
)

# Condition 4: Equal ERS, different theta
cat("  4. equal_ers:\n")
cat("     - Group 1: theta~N(0,1), ERS~N(0.5,1)\n")
cat("     - Group 2: theta~N(0.3,1), ERS~N(0.5,1)\n")
cat("     - Purpose: ERS equal, theta different (paradox scenario)\n")
cat("\n")

equal_ers <- data.frame(
  theta_mean = c(0, 0.3),
  theta_sd = c(1, 1),
  ers_mean = c(0.5, 0.5),
  ers_sd = c(1, 1)
)

# Combine into list
ers_conditions <- list(
  low_ers = low_ers,
  medium_ers = medium_ers,
  high_ers_diff = high_ers_diff,
  equal_ers = equal_ers
)


# 4. Generate All Conditions ----

cat("=======================================================\n")
cat("  Starting Data Generation                            \n")
cat("=======================================================\n")
cat("\n")

all_datasets <- list()
start_time_total <- Sys.time()

for (ers_name in names(ers_conditions)) {
  for (n_size in SAMPLE_SIZES) {

    condition_name <- paste0(ers_name, "_N", n_size)

    cat("\n")
    cat("--- Condition:", condition_name, "---\n")
    cat("  ERS parameters:", ers_name, "\n")
    cat("  Sample size per group:", n_size, "\n")
    cat("  Generating", N_REPLICATIONS, "replications...\n")

    start_time_condition <- Sys.time()

    # Generate replications
    replications <- vector("list", N_REPLICATIONS)

    for (rep in 1:N_REPLICATIONS) {

      # Progress indicator
      if (rep %% 100 == 0) {
        cat("    Replication", rep, "/", N_REPLICATIONS, "\n")
      }

      # Generate data
      replications[[rep]] <- generate_mnrm_falk(
        N = n_size,
        n_items = N_ITEMS,
        n_groups = 2,
        categories = CATEGORIES,
        alpha = ALPHA,
        thresholds = THRESHOLDS_AVERAGE,
        theta_df = ers_conditions[[ers_name]]
      )
    }

    # Store in main list
    all_datasets[[condition_name]] <- replications

    # Report timing
    time_elapsed <- difftime(Sys.time(), start_time_condition, units = "secs")
    cat("  ✓ Condition completed in", round(time_elapsed, 1), "seconds\n")
  }
}


# 5. Save Results ----

cat("\n")
cat("=======================================================\n")
cat("  Saving Results                                      \n")
cat("=======================================================\n")
cat("\n")

# Save all datasets
saveRDS(all_datasets, file = "data/simulated/all_conditions.rds")
cat("  ✓ All datasets saved to: data/simulated/all_conditions.rds\n")

# Create summary table
summary_table <- data.frame(
  condition = names(all_datasets),
  n_replications = sapply(all_datasets, length),
  sample_size_per_group = rep(SAMPLE_SIZES, each = 4),
  total_sample_size = rep(SAMPLE_SIZES * 2, each = 4)
)

write.csv(summary_table, file = "output/tables/simulation_summary.csv", row.names = FALSE)
cat("  ✓ Summary table saved to: output/tables/simulation_summary.csv\n")


# 6. Verification ----

cat("\n")
cat("=======================================================\n")
cat("  Verification                                        \n")
cat("=======================================================\n")
cat("\n")

# Check one example dataset
example_condition <- "medium_ers_N500"
example_rep <- all_datasets[[example_condition]][[1]]

cat("Example dataset (", example_condition, ", replication 1):\n", sep = "")
cat("\n")
cat("Data dimensions:", nrow(example_rep$data), "rows ×", ncol(example_rep$data), "columns\n")
cat("\n")

cat("First 6 rows:\n")
print(head(example_rep$data))
cat("\n")

cat("True parameter means by group:\n")
cat("\n")

group_summary <- example_rep$data %>%
  group_by(country) %>%
  summarise(
    n = n(),
    mean_true_theta = mean(true_theta),
    sd_true_theta = sd(true_theta),
    mean_true_ers = mean(true_ers),
    sd_true_ers = sd(true_ers),
    .groups = "drop"
  )

print(group_summary)
cat("\n")

# Calculate observed score means
observed_summary <- example_rep$data %>%
  select(country, starts_with("Item")) %>%
  pivot_longer(cols = starts_with("Item"), names_to = "item", values_to = "response") %>%
  group_by(country) %>%
  summarise(
    mean_response = mean(response, na.rm = TRUE),
    sd_response = sd(response, na.rm = TRUE),
    .groups = "drop"
  )

cat("Observed response means by group:\n")
print(observed_summary)
cat("\n")


# 7. Final Summary ----

total_time <- difftime(Sys.time(), start_time_total, units = "mins")

cat("=======================================================\n")
cat("  Data Generation Completed Successfully!             \n")
cat("=======================================================\n")
cat("\n")
cat("Total time:", round(total_time, 2), "minutes\n")
cat("Total datasets generated:", sum(summary_table$n_replications), "\n")
cat("Total participants:", sum(summary_table$n_replications * summary_table$total_sample_size), "\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Run: source('R/simulation/02_correction_methods.R')\n")
cat("  2. Or run: source('R/run_all_analyses.R')\n")
cat("\n")
cat("=======================================================\n")
