# ============================================================================
# Script: 05_debug_na_values.R
# Purpose: Debug NA values in Z-score and IRTree results
# Author: ERS Research Project
# Date: 2025-11-30
# ============================================================================

library(tidyverse)

cat("\n")
cat("=======================================================\n")
cat("  DEBUGGING NA VALUES IN CORRECTION METHODS           \n")
cat("=======================================================\n")
cat("\n")

# Load data
cat("Loading simulation data...\n")
all_datasets <- readRDS("data/simulated/all_conditions_full_factorial.rds")

# Test with first condition
condition_name <- names(all_datasets)[1]
cat("Testing with condition:", condition_name, "\n\n")

# Get first replication
rep_data <- all_datasets[[condition_name]][[1]]
sim_data <- rep_data$data

# Get item columns
item_cols <- grep("^Item", names(sim_data), value = TRUE)
n_items <- length(item_cols)

cat("Number of items:", n_items, "\n")
cat("Number of rows:", nrow(sim_data), "\n")
cat("\n")

# Source functions
source("R/functions/correction_methods.R")
source("R/functions/ers_indices.R")

# Test 1: Z-score correction
cat("=== Testing Z-score Correction ===\n")
cat("\n")

# Check for zero variance cases
item_data <- sim_data[, item_cols]
person_sd <- apply(item_data, 1, sd, na.rm = TRUE)

cat("Person SD statistics:\n")
cat("  Min:", min(person_sd, na.rm = TRUE), "\n")
cat("  Max:", max(person_sd, na.rm = TRUE), "\n")
cat("  N with SD=0:", sum(person_sd == 0, na.rm = TRUE), "\n")
cat("  N with NA SD:", sum(is.na(person_sd)), "\n")
cat("\n")

# Apply z-score
data_zscore <- apply_zscore_correction(sim_data, items = item_cols, suffix = "_z")
zscore_cols <- paste0(item_cols, "_z")

# Check for NAs in z-scored data
cat("Z-scored data check:\n")
cat("  Columns created:", length(zscore_cols), "\n")
cat("  NA counts in z-scored columns:\n")

for (col in zscore_cols[1:min(5, length(zscore_cols))]) {
  n_na <- sum(is.na(data_zscore[[col]]))
  cat("    ", col, ":", n_na, "NAs\n")
}
cat("\n")

# Calculate mean z-score
zscore_means <- rowMeans(data_zscore[, zscore_cols], na.rm = TRUE)
cat("Mean z-score statistics:\n")
cat("  N total:", length(zscore_means), "\n")
cat("  N NA:", sum(is.na(zscore_means)), "\n")
cat("  N Inf:", sum(is.infinite(zscore_means)), "\n")
cat("  N NaN:", sum(is.nan(zscore_means)), "\n")
cat("\n")

# Group-level means
cat("Group-level mean z-scores:\n")
group_means_zscore <- data_zscore %>%
  mutate(zscore_mean = rowMeans(across(all_of(zscore_cols)), na.rm = TRUE)) %>%
  group_by(country) %>%
  summarise(
    mean_zscore = mean(zscore_mean, na.rm = TRUE),
    n = n(),
    n_na = sum(is.na(zscore_mean))
  )
print(group_means_zscore)
cat("\n")


# Test 2: Covariate method
cat("=== Testing Covariate Method ===\n")
cat("\n")

# Add ERS index
data_ers <- add_ers_index(sim_data, items = item_cols, index_name = "ers_index")

cat("ERS index statistics:\n")
cat("  N total:", nrow(data_ers), "\n")
cat("  N NA:", sum(is.na(data_ers$ers_index)), "\n")
cat("  Min:", min(data_ers$ers_index, na.rm = TRUE), "\n")
cat("  Max:", max(data_ers$ers_index, na.rm = TRUE), "\n")
cat("\n")

# Calculate raw scale score
data_ers$scale_score_raw <- rowMeans(data_ers[, item_cols], na.rm = TRUE)

# Apply covariate correction
data_covariate <- apply_covariate_correction(
  data = data_ers,
  items = item_cols,
  ers_index = "ers_index",
  score_method = "mean",
  correction_method = "residuals"
)

cat("Covariate-corrected scores:\n")
if ("scale_score_ers_corrected" %in% names(data_covariate)) {
  cat("  N total:", length(data_covariate$scale_score_ers_corrected), "\n")
  cat("  N NA:", sum(is.na(data_covariate$scale_score_ers_corrected)), "\n")
} else {
  cat("  ERROR: scale_score_ers_corrected column not found!\n")
  cat("  Available columns:", paste(names(data_covariate), collapse = ", "), "\n")
}
cat("\n")


# Test 3: Check raw data quality
cat("=== Checking Raw Data Quality ===\n")
cat("\n")

# Check for missing values in raw data
cat("Missing values in raw item data:\n")
for (col in item_cols[1:min(5, length(item_cols))]) {
  n_na <- sum(is.na(sim_data[[col]]))
  cat("  ", col, ":", n_na, "NAs\n")
}
cat("\n")

# Check response distributions
cat("Response value distributions (first 3 items):\n")
for (col in item_cols[1:min(3, length(item_cols))]) {
  cat("  ", col, ":\n")
  print(table(sim_data[[col]], useNA = "ifany"))
  cat("\n")
}


# Test 4: Sample diagnostic cases
cat("=== Sample Diagnostic Cases ===\n")
cat("\n")

# Find cases with all same responses (potential issue)
same_response_cases <- which(person_sd == 0)

if (length(same_response_cases) > 0) {
  cat("Cases with all same responses (SD=0):\n")
  cat("  N cases:", length(same_response_cases), "\n")
  cat("\n")
  cat("  Example case (first 5 items):\n")
  print(sim_data[same_response_cases[1], item_cols[1:min(5, length(item_cols))]])
  cat("\n")
  cat("  Z-scores for this case:\n")
  print(data_zscore[same_response_cases[1], zscore_cols[1:min(5, length(zscore_cols))]])
  cat("\n")
} else {
  cat("No cases with SD=0 found (good!)\n")
  cat("\n")
}


cat("=======================================================\n")
cat("  Debug Analysis Complete                             \n")
cat("=======================================================\n")
cat("\n")
cat("Summary:\n")
cat("  1. Check if z-score NAs are due to SD=0 cases\n")
cat("  2. Verify that correction functions are working\n")
cat("  3. Investigate any unexpected NA patterns\n")
cat("\n")
cat("Next steps:\n")
cat("  - If NAs persist, check correction_methods.R functions\n")
cat("  - Verify data generation is producing valid responses\n")
cat("  - Check for numeric precision issues\n")
cat("\n")
