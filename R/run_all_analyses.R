# ============================================================================
# Script: run_all_analyses.R
# Purpose: Master script to run entire ERS analysis pipeline
# Author: ERS Research Project
# Date: 2025-11-28
# ============================================================================
#
# This script runs the complete ERS analysis in the correct order:
# 1. Simulation study (data generation, correction, evaluation)
# 2. PISA data analysis (preparation, models, visualization)
#
# Usage:
#   source("R/run_all_analyses.R")
#
# Or run individual phases:
#   RUN_SIMULATION <- TRUE
#   RUN_PISA <- FALSE
#   source("R/run_all_analyses.R")
# ============================================================================

# Configuration ----

# Set which analyses to run
if (!exists("RUN_SIMULATION")) RUN_SIMULATION <- TRUE
if (!exists("RUN_PISA")) RUN_PISA <- TRUE
if (!exists("RUN_VISUALIZATIONS")) RUN_VISUALIZATIONS <- TRUE

cat("\n")
cat("=======================================================\n")
cat("  ERS CORRECTION METHODS ANALYSIS - MASTER SCRIPT     \n")
cat("=======================================================\n")
cat("\n")
cat("Configuration:\n")
cat("  Run Simulation Study:    ", RUN_SIMULATION, "\n")
cat("  Run PISA Analysis:       ", RUN_PISA, "\n")
cat("  Run Visualizations:      ", RUN_VISUALIZATIONS, "\n")
cat("\n")


# Create output directories ----

required_dirs <- c(
  "data/raw",
  "data/processed",
  "data/simulated",
  "output/figures",
  "output/tables",
  "output/reports",
  "output/models"
)

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

cat("\n")


# ============================================================================
# PART 1: SIMULATION STUDY ----
# ============================================================================

if (RUN_SIMULATION) {

  cat("=======================================================\n")
  cat("PART 1: SIMULATION STUDY                              \n")
  cat("=======================================================\n\n")

  # 1.1 Data Generation ----
  cat("--- Step 1: Generating simulated data ---\n")
  start_time <- Sys.time()

  tryCatch({
    source("R/simulation/01_data_generation.R")
    cat("✓ Data generation completed\n")
  }, error = function(e) {
    cat("✗ Error in data generation:", e$message, "\n")
  })

  cat("Time elapsed:", difftime(Sys.time(), start_time, units = "secs"), "seconds\n\n")


  # 1.2 Correction Methods ----
  cat("--- Step 2: Applying ERS correction methods ---\n")
  start_time <- Sys.time()

  tryCatch({
    source("R/simulation/02_correction_methods.R")
    cat("✓ Correction methods applied\n")
  }, error = function(e) {
    cat("✗ Error in correction methods:", e$message, "\n")
  })

  cat("Time elapsed:", difftime(Sys.time(), start_time, units = "secs"), "seconds\n\n")


  # 1.3 Evaluation ----
  cat("--- Step 3: Evaluating methods (Bias, RMSE) ---\n")
  start_time <- Sys.time()

  tryCatch({
    source("R/simulation/03_evaluation.R")
    cat("✓ Evaluation completed\n")
  }, error = function(e) {
    cat("✗ Error in evaluation:", e$message, "\n")
  })

  cat("Time elapsed:", difftime(Sys.time(), start_time, units = "secs"), "seconds\n\n")

  cat("PART 1 COMPLETED: Simulation study\n\n")
}


# ============================================================================
# PART 2: PISA DATA ANALYSIS ----
# ============================================================================

if (RUN_PISA) {

  cat("=======================================================\n")
  cat("PART 2: PISA DATA ANALYSIS                            \n")
  cat("=======================================================\n\n")

  # 2.1 Data Preparation ----
  cat("--- Step 1: Preparing PISA data ---\n")
  start_time <- Sys.time()

  tryCatch({
    source("R/analysis/01_data_preparation.R")
    cat("✓ Data preparation completed\n")
  }, error = function(e) {
    cat("✗ Error in data preparation:", e$message, "\n")
    cat("  Note: If PISA data file is not available, placeholder data will be used.\n")
  })

  cat("Time elapsed:", difftime(Sys.time(), start_time, units = "secs"), "seconds\n\n")


  # 2.2 ERS Measurement and Models ----
  cat("--- Step 2: Calculating ERS indices and fitting all models ---\n")
  start_time <- Sys.time()

  tryCatch({
    source("R/analysis/02_ers_measurement_and_models.R")
    cat("✓ ERS measurement and models completed\n")
  }, error = function(e) {
    cat("✗ Error in ERS measurement/models:", e$message, "\n")
  })

  cat("Time elapsed:", difftime(Sys.time(), start_time, units = "mins"), "minutes\n\n")

  cat("PART 2 COMPLETED: PISA analysis\n\n")
}


# ============================================================================
# PART 3: VISUALIZATIONS ----
# ============================================================================

if (RUN_VISUALIZATIONS) {

  cat("=======================================================\n")
  cat("PART 3: CREATING VISUALIZATIONS                       \n")
  cat("=======================================================\n\n")

  # 3.1 Country Comparisons ----
  if (RUN_PISA) {
    cat("--- Creating PISA country comparison plots ---\n")
    start_time <- Sys.time()

    tryCatch({
      source("R/visualization/country_comparisons.R")
      cat("✓ PISA visualizations created\n")
    }, error = function(e) {
      cat("✗ Error in PISA visualizations:", e$message, "\n")
    })

    cat("Time elapsed:", difftime(Sys.time(), start_time, units = "secs"),
       "seconds\n\n")
  }

  cat("PART 3 COMPLETED: Visualizations\n\n")
}


# ============================================================================
# FINAL SUMMARY ----
# ============================================================================

cat("\n")
cat("=======================================================\n")
cat("  ALL ANALYSES COMPLETED!                              \n")
cat("=======================================================\n\n")

cat("Output locations:\n")
cat("  - Simulated data:        data/simulated/\n")
cat("  - Processed PISA data:   data/processed/\n")
cat("  - Tables:                output/tables/\n")
cat("  - Figures:               output/figures/\n")
cat("  - Model objects:         output/models/\n")

cat("\n")
cat("Key output files:\n")

if (RUN_SIMULATION) {
  cat("\nSimulation Study:\n")
  cat("  - simulation_evaluation_metrics.csv\n")
  cat("  - simulation_ranking_changes.csv\n")
  cat("  - simulation_rmse_comparison.png\n")
  cat("  - simulation_bias_comparison.png\n")
}

if (RUN_PISA) {
  cat("\nPISA Analysis:\n")
  cat("  - pisa_all_models_results.csv\n")
  cat("  - country_means_comparison.png\n")
  cat("  - ranking_changes.png\n")
  cat("  - ers_vs_country_means.png\n")
  cat("  - rank_changes_heatmap.png\n")
  cat("  - method_agreement.png\n")
}

cat("\n")
cat("=======================================================\n")
cat("Analysis pipeline completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=======================================================\n\n")

# Save session info ----
sink("output/reports/session_info.txt")
cat("ERS Analysis - Session Information\n")
cat("===================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
sink()

cat("Session info saved to: output/reports/session_info.txt\n\n")
