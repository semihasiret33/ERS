# ============================================================================
# Script: 00_irtree_calibration.R
# Purpose: Calibrate IRTree parameters from MNRM data
# Author: ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# This script calibrates IRTree item parameters by:
# 1. Generating a large MNRM dataset (N=50,000 per group)
# 2. Fitting an IRTree model to this data
# 3. Extracting item parameters to use as "true" IRTree parameters
#
# This approach ensures maximum comparability between MNRM and IRTree models
# following Schoenmakers et al. (2023) methodology.
# ============================================================================

library(mirt)
library(tidyverse)

# Source the MNRM generator and IRTree functions
source("R/simulation/00_mnrm_falk_generator.R")
source("R/functions/irtree_functions.R")


#' Calibrate IRTree Parameters from MNRM Data
#'
#' Generates large MNRM dataset and fits IRTree model to obtain
#' IRTree item parameters for use in data generation.
#'
#' @param N_calibration Sample size per group for calibration (default = 50000)
#' @param n_items Number of items
#' @param alpha MNRM discrimination parameters
#' @param thresholds MNRM thresholds
#' @param theta_df Group parameters for calibration data
#' @return List containing IRTree node parameters
calibrate_irtree_from_mnrm <- function(N_calibration = 50000,
                                       n_items = 10,
                                       alpha = c(1.5, 1.5),
                                       thresholds = c(-1, 0, 1),
                                       theta_df) {

  cat("\n===============================================\n")
  cat("  IRTree Parameter Calibration from MNRM      \n")
  cat("===============================================\n\n")

  # Step 1: Generate large MNRM dataset
  cat("Step 1: Generating large MNRM dataset (N =", N_calibration, "per group)...\n")

  calibration_data <- generate_mnrm_falk(
    N = N_calibration,
    n_items = n_items,
    alpha = alpha,
    thresholds = thresholds,
    theta_df = theta_df
  )

  cat("  MNRM data generated successfully!\n")

  # Step 2: Convert to IRTree pseudo-item format
  cat("\nStep 2: Converting to IRTree pseudo-item format...\n")

  item_cols <- paste0("Item", 1:n_items)
  obs_matrix <- as.matrix(calibration_data$data[, item_cols])
  n_total <- nrow(obs_matrix)

  # Create pseudo-item matrix (3 nodes per item)
  tree_data <- matrix(NA, nrow = n_total, ncol = n_items * 3)

  for (i in 1:n_items) {
    # Node 1: Agree (1) vs Disagree (0)
    tree_data[obs_matrix[, i] <= 2, (i - 1) * 3 + 1] <- 0  # Disagree
    tree_data[obs_matrix[, i] >= 3, (i - 1) * 3 + 1] <- 1  # Agree

    # Node 2: Given disagree, extreme (0) vs moderate (1)
    tree_data[obs_matrix[, i] == 2, (i - 1) * 3 + 2] <- 1  # Moderate
    tree_data[obs_matrix[, i] == 1, (i - 1) * 3 + 2] <- 0  # Extreme

    # Node 3: Given agree, moderate (0) vs extreme (1)
    tree_data[obs_matrix[, i] == 3, (i - 1) * 3 + 3] <- 0  # Moderate
    tree_data[obs_matrix[, i] == 4, (i - 1) * 3 + 3] <- 1  # Extreme
  }

  # Name the pseudo-items
  items <- sort(rep(1:n_items, 3))
  items <- paste0("item", items)
  nodes <- 1:3
  nodes <- rep(paste0("node", nodes), n_items)
  colnames(tree_data) <- paste0(items, nodes)

  cat("  Pseudo-items created:", ncol(tree_data), "nodes\n")

  # Step 3: Fit IRTree model
  cat("\nStep 3: Fitting IRTree model to calibration data...\n")
  cat("  This may take several minutes...\n")

  # Initialize custom IRTree nodes
  irtree_nodes <- initialize_irtree_nodes()

  # Create group variable
  group_var <- calibration_data$data$country

  # Determine which nodes load on ERS (nodes 2 and 3)
  ers_loadings <- rep(NA, n_items * 2)
  for (i in 1:n_items) {
    ers_loadings[(i - 1) * 2 + 1] <- (i - 1) * 3 + 2
    ers_loadings[(i - 1) * 2 + 2] <- (i - 1) * 3 + 3
  }

  # Specify item types
  itemtypes <- rep(c("2PL", "Node2TreePL", "Node3TreePL"), n_items)

  # Build mirt model syntax
  tree_n_items <- ncol(tree_data)

  model_syntax <- paste0(
    "Theta = 1-", tree_n_items, "\n",
    "ERS = ", paste0(ers_loadings, collapse = ","), "\n",
    "CONSTRAINB = (1-", tree_n_items, ", a1), (1-", tree_n_items, ", a2), ",
    "(1-", tree_n_items, ", d)\n",
    "FREE = (GROUP, COV_21)"
  )

  # Fit the model
  tryCatch({
    irtree_fit <- multipleGroup(
      data = tree_data,
      model = model_syntax,
      group = factor(group_var),
      itemtype = itemtypes,
      method = "EM",
      invariance = c("free_mean", "free_var"),
      technical = list(NCYCLES = 2000),
      verbose = TRUE,
      customItems = list(
        Node2TreePL = irtree_nodes$Node2TreePL,
        Node3TreePL = irtree_nodes$Node3TreePL
      )
    )

    cat("\n  IRTree model fitted successfully!\n")

    # Step 4: Extract parameters
    cat("\nStep 4: Extracting IRTree parameters...\n")

    coef_list <- coef(irtree_fit, simplify = TRUE)

    # Extract node parameters for each item
    node1_pars <- matrix(NA, nrow = n_items, ncol = 3)
    node2_pars <- matrix(NA, nrow = n_items, ncol = 3)
    node3_pars <- matrix(NA, nrow = n_items, ncol = 3)

    for (i in 1:n_items) {
      # Node 1: 2PL parameters (a1, a2, d)
      node1_pars[i, ] <- coef_list$Group_1[[(i - 1) * 3 + 1]][1:3]

      # Node 2: Custom node parameters
      node2_pars[i, ] <- coef_list$Group_1[[(i - 1) * 3 + 2]]

      # Node 3: Custom node parameters
      node3_pars[i, ] <- coef_list$Group_1[[(i - 1) * 3 + 3]]
    }

    colnames(node1_pars) <- c("a1", "a2", "d")
    colnames(node2_pars) <- c("a1", "a2", "d")
    colnames(node3_pars) <- c("a1", "a2", "d")

    cat("\n  Parameters extracted successfully!\n")

    # Print parameter summaries
    cat("\nNode 1 parameters (Agree vs Disagree):\n")
    print(summary(node1_pars))

    cat("\nNode 2 parameters (Extreme vs Moderate Disagree):\n")
    print(summary(node2_pars))

    cat("\nNode 3 parameters (Moderate vs Extreme Agree):\n")
    print(summary(node3_pars))

    # Create result object
    result <- list(
      node1_pars = node1_pars,
      node2_pars = node2_pars,
      node3_pars = node3_pars,
      irtree_fit = irtree_fit,
      calibration_data = calibration_data,
      n_items = n_items,
      alpha_mnrm = alpha,
      thresholds_mnrm = thresholds
    )

    cat("\n===============================================\n")
    cat("  IRTree Calibration Completed Successfully!  \n")
    cat("===============================================\n\n")

    return(result)

  }, error = function(e) {
    cat("\nERROR: IRTree fitting failed!\n")
    cat("Error message:", e$message, "\n")
    cat("\nReturning NULL. You may need to:\n")
    cat("  1. Reduce calibration sample size\n")
    cat("  2. Simplify model constraints\n")
    cat("  3. Check data quality\n")

    return(NULL)
  })
}


#' Calibrate IRTree for Multiple Conditions
#'
#' Calibrates IRTree parameters for both threshold sets
#' (average and difficult items)
#'
#' @param save_results Whether to save results to RDS file
#' @return List with calibration results for both conditions
calibrate_all_irtree_parameters <- function(save_results = TRUE) {

  # Define calibration data parameters
  # Use neutral ERS (mean = 0 for both groups) for calibration
  theta_df_calibration <- data.frame(
    theta_mean = c(0, 0),
    theta_sd = c(1, 1),
    ers_mean = c(0, 0),
    ers_sd = c(1, 1)
  )

  # Calibrate for average thresholds (10 items)
  cat("\n\n")
  cat("====================================================\n")
  cat("  Calibrating IRTree: 10 items, Average difficulty \n")
  cat("====================================================\n")

  irtree_params_10_avg <- calibrate_irtree_from_mnrm(
    N_calibration = 50000,
    n_items = 10,
    alpha = c(1.5, 1.5),
    thresholds = c(-1, 0, 1),
    theta_df = theta_df_calibration
  )

  # Calibrate for difficult thresholds (10 items)
  cat("\n\n")
  cat("====================================================\n")
  cat("  Calibrating IRTree: 10 items, Difficult items    \n")
  cat("====================================================\n")

  irtree_params_10_diff <- calibrate_irtree_from_mnrm(
    N_calibration = 50000,
    n_items = 10,
    alpha = c(1.5, 1.5),
    thresholds = c(0, 1, 2),
    theta_df = theta_df_calibration
  )

  # Combine results
  all_params <- list(
    items_10_average = irtree_params_10_avg,
    items_10_difficult = irtree_params_10_diff
  )

  # Save if requested
  if (save_results) {
    saveRDS(all_params, file = "data/simulated/irtree_calibration_parameters.rds")
    cat("\nCalibration parameters saved to: data/simulated/irtree_calibration_parameters.rds\n")
  }

  return(all_params)
}


# Example usage (uncomment to run):
# calibration_results <- calibrate_all_irtree_parameters(save_results = TRUE)
