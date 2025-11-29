# ============================================================================
# Script: 01_data_generation.R
# Purpose: Generate simulation data for all conditions (Full Factorial Design)
# Author: ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# FULL FACTORIAL SIMULATION DESIGN:
# - 4 ERS conditions (low, medium, high_diff, equal)
# - 3 Sample sizes (N=100, N=500, N=1000)
# - 3 Scale lengths (10, 20, 30 items)
# - 2 Likert types (4-point, 5-point)
# - 100 replications per condition
#
# Total: 4 × 3 × 3 × 2 = 72 conditions × 100 replications = 7,200 datasets
#
# Data generation uses MNRM-Falk approach (Falk & Cai, 2016) with:
# - Alpha = 1.5 for both dimensions
# - Varying item difficulties (m values)
# - Thresholds adjusted for scale length
# ============================================================================

# 1. Setup ----
library(MASS)
library(tidyverse)

# Source data generation functions
source("R/simulation/00_mnrm_falk_generator.R")

set.seed(12345)

cat("\n")
cat("=======================================================\n")
cat("  ERS FULL FACTORIAL SIMULATION                       \n")
cat("=======================================================\n")
cat("\n")


# 2. Define Simulation Parameters ----

cat("Full Factorial Design Parameters:\n")
cat("  - ERS conditions: 4\n")
cat("  - Sample sizes: 3 (100, 500, 1000)\n")
cat("  - Scale lengths: 3 (10, 20, 30 items)\n")
cat("  - Likert types: 2 (4-point, 5-point)\n")
cat("  - Replications per condition: 100\n")
cat("\n")
cat("  Total conditions: 72\n")
cat("  Total datasets: 7,200\n")
cat("\n")

# Fixed parameters
ALPHA <- c(1.5, 1.5)
N_REPLICATIONS <- 100

# Factorial design factors
SAMPLE_SIZES <- c(100, 500, 1000)
SCALE_LENGTHS <- c(10, 20, 30)
LIKERT_TYPES <- c(4, 5)

# Thresholds for different Likert scales
THRESHOLDS_4POINT <- c(-1, 0, 1)
THRESHOLDS_5POINT <- c(-1.5, -0.5, 0.5, 1.5)


# 3. Define ERS Conditions ----

cat("ERS Conditions:\n")
cat("\n")

# Condition 1: Low ERS (both groups low, no difference)
cat("  1. low_ers: Both groups ERS~N(0,1)\n")
low_ers <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 0),
  ers_sd = c(1, 1)
)

# Condition 2: Medium ERS (group difference of 0.5 SD)
cat("  2. medium_ers: Group 2 ERS~N(0.5,1)\n")
medium_ers <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 0.5),
  ers_sd = c(1, 1)
)

# Condition 3: High ERS difference (1 SD difference)
cat("  3. high_ers_diff: Group 2 ERS~N(1.0,1)\n")
high_ers_diff <- data.frame(
  theta_mean = c(0, 0),
  theta_sd = c(1, 1),
  ers_mean = c(0, 1.0),
  ers_sd = c(1, 1)
)

# Condition 4: Equal ERS, different theta
cat("  4. equal_ers: Both ERS~N(0.5,1), Group 2 theta~N(0.3,1)\n")
equal_ers <- data.frame(
  theta_mean = c(0, 0.3),
  theta_sd = c(1, 1),
  ers_mean = c(0.5, 0.5),
  ers_sd = c(1, 1)
)

cat("\n")

# Combine into list
ers_conditions <- list(
  low_ers = low_ers,
  medium_ers = medium_ers,
  high_ers_diff = high_ers_diff,
  equal_ers = equal_ers
)


# 4. Generate All Conditions ----

cat("=======================================================\n")
cat("  Starting Full Factorial Data Generation             \n")
cat("=======================================================\n")
cat("\n")

all_datasets <- list()
start_time_total <- Sys.time()
condition_counter <- 0

# Nested loops for full factorial design
for (ers_name in names(ers_conditions)) {
  for (n_size in SAMPLE_SIZES) {
    for (n_items in SCALE_LENGTHS) {
      for (n_categories in LIKERT_TYPES) {

        condition_counter <- condition_counter + 1

        # Create condition name
        condition_name <- paste0(
          ers_name,
          "_N", n_size,
          "_I", n_items,
          "_L", n_categories
        )

        cat("\n")
        cat("--- Condition", condition_counter, "/72:", condition_name, "---\n")
        cat("  ERS:", ers_name, "\n")
        cat("  Sample size:", n_size, "per group\n")
        cat("  Items:", n_items, "\n")
        cat("  Likert:", n_categories, "-point\n")
        cat("  Generating", N_REPLICATIONS, "replications...\n")

        start_time_condition <- Sys.time()

        # Select appropriate thresholds
        if (n_categories == 4) {
          thresholds_use <- THRESHOLDS_4POINT
        } else {
          thresholds_use <- THRESHOLDS_5POINT
        }

        # Generate replications
        replications <- vector("list", N_REPLICATIONS)

        for (rep in 1:N_REPLICATIONS) {

          # Progress indicator (every 25 reps)
          if (rep %% 25 == 0) {
            cat("    Rep", rep, "/", N_REPLICATIONS, "\n")
          }

          # Generate data
          replications[[rep]] <- generate_mnrm_falk(
            N = n_size,
            n_items = n_items,
            n_groups = 2,
            categories = n_categories,
            alpha = ALPHA,
            thresholds = thresholds_use,
            theta_df = ers_conditions[[ers_name]]
          )
        }

        # Store in main list
        all_datasets[[condition_name]] <- replications

        # Report timing
        time_elapsed <- difftime(Sys.time(), start_time_condition, units = "secs")
        cat("  ✓ Completed in", round(time_elapsed, 1), "seconds\n")
      }
    }
  }
}


# 5. Save Results ----

cat("\n")
cat("=======================================================\n")
cat("  Saving Results                                      \n")
cat("=======================================================\n")
cat("\n")

# Save all datasets
saveRDS(all_datasets, file = "data/simulated/all_conditions_full_factorial.rds")
cat("  ✓ All datasets saved to: data/simulated/all_conditions_full_factorial.rds\n")

# Create detailed summary table
summary_table <- data.frame(
  condition = names(all_datasets),
  n_replications = sapply(all_datasets, length)
)

# Parse condition names to extract factors
summary_table <- summary_table %>%
  mutate(
    ers_condition = str_extract(condition, "^[^_]+"),
    sample_size = as.numeric(str_extract(condition, "(?<=_N)\\d+")),
    n_items = as.numeric(str_extract(condition, "(?<=_I)\\d+")),
    likert_type = as.numeric(str_extract(condition, "(?<=_L)\\d+")),
    total_n = sample_size * 2
  )

write.csv(summary_table, file = "output/tables/simulation_summary_full_factorial.csv",
          row.names = FALSE)
cat("  ✓ Summary table saved to: output/tables/simulation_summary_full_factorial.csv\n")


# 6. Verification ----

cat("\n")
cat("=======================================================\n")
cat("  Verification                                        \n")
cat("=======================================================\n")
cat("\n")

# Check an example dataset from each factorial level
example_conditions <- c(
  "medium_ers_N500_I10_L4",   # Baseline
  "medium_ers_N500_I20_L4",   # Long scale
  "medium_ers_N500_I10_L5"    # 5-point Likert
)

for (cond in example_conditions) {
  if (cond %in% names(all_datasets)) {
    example_rep <- all_datasets[[cond]][[1]]

    cat("\n", cond, ":\n", sep = "")
    cat("  Dimensions:", nrow(example_rep$data), "rows ×",
        ncol(example_rep$data), "columns\n")

    group_summary <- example_rep$data %>%
      group_by(country) %>%
      summarise(
        n = n(),
        mean_theta = mean(true_theta),
        mean_ers = mean(true_ers),
        .groups = "drop"
      )

    print(group_summary)
  }
}


# 7. Create Factorial Summary Statistics ----

cat("\n")
cat("=======================================================\n")
cat("  Factorial Design Summary                            \n")
cat("=======================================================\n")
cat("\n")

factorial_summary <- summary_table %>%
  group_by(ers_condition, sample_size, n_items, likert_type) %>%
  summarise(
    n_conditions = n(),
    total_datasets = sum(n_replications),
    .groups = "drop"
  )

cat("Conditions by ERS type:\n")
print(table(summary_table$ers_condition))
cat("\n")

cat("Conditions by sample size:\n")
print(table(summary_table$sample_size))
cat("\n")

cat("Conditions by number of items:\n")
print(table(summary_table$n_items))
cat("\n")

cat("Conditions by Likert type:\n")
print(table(summary_table$likert_type))
cat("\n")


# 8. Final Summary ----

total_time <- difftime(Sys.time(), start_time_total, units = "mins")

cat("=======================================================\n")
cat("  Full Factorial Data Generation Completed!           \n")
cat("=======================================================\n")
cat("\n")
cat("Total conditions:", nrow(summary_table), "\n")
cat("Total replications:", sum(summary_table$n_replications), "\n")
cat("Total datasets:", sum(summary_table$n_replications), "\n")
cat("Total time:", round(total_time, 2), "minutes\n")
cat("Average time per condition:", round(total_time / 72, 2), "minutes\n")
cat("\n")

# Estimate file sizes
datasets_generated <- sum(summary_table$n_replications)
total_participants <- sum(summary_table$n_replications * summary_table$total_n)
cat("Total participants simulated:", format(total_participants, big.mark = ","), "\n")
cat("\n")

cat("Next steps:\n")
cat("  1. Run: source('R/simulation/02_correction_methods.R')\n")
cat("  2. Or run: source('R/run_all_analyses.R')\n")
cat("\n")
cat("=======================================================\n")
