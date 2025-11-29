# ============================================================================
# Script: 00_mnrm_falk_generator.R
# Purpose: MNRM-Falk data generation function (based on Falk & Cai 2016)
# Author: ERS Research Project
# Date: 2025-11-29
# ============================================================================
#
# This script implements the MNRM data generation approach from:
# - Falk & Cai (2016) - Multidimensional Nominal Response Model
# - Schoenmakers et al. (2023) - Model comparison study
#
# Key features:
# - Alpha = 1.5 for both theta and ERS dimensions
# - Varying item difficulties using m values
# - Two threshold sets (average vs difficult)
# - ERS contamination through s-matrix scoring
# ============================================================================

library(MASS)  # For mvrnorm


#' Generate MNRM-Falk Data
#'
#' Generates multigroup data under the Multidimensional Nominal Response Model
#' following Falk & Cai (2016) methodology.
#'
#' @param N Sample size per group
#' @param n_items Number of items
#' @param n_groups Number of groups (default = 2 for cross-cultural comparison)
#' @param categories Number of response categories (default = 4)
#' @param alpha Item discrimination parameters (length 2: theta, ERS)
#' @param thresholds Base item thresholds (length 3 for 4 categories)
#' @param theta_df Data frame with group parameters (theta_mean, theta_sd, ers_mean, ers_sd)
#' @return List containing response data, true parameters, group membership
#' @references Falk & Cai (2016), Schoenmakers et al. (2023)
generate_mnrm_falk <- function(N = 500,
                                n_items = 10,
                                n_groups = 2,
                                categories = 4,
                                alpha = c(1.5, 1.5),
                                thresholds = c(-1, 0, 1),
                                theta_df) {

  # Validate inputs
  if (nrow(theta_df) != n_groups) {
    stop("theta_df must have one row per group")
  }

  if (!all(c("theta_mean", "theta_sd", "ers_mean", "ers_sd") %in% names(theta_df))) {
    stop("theta_df must contain: theta_mean, theta_sd, ers_mean, ers_sd")
  }

  # Generate s-matrix for MNRM
  # Content dimension: linear scoring 0, 1, 2, 3
  s_content <- matrix(0:(categories - 1), nrow = 1)

  # ERS dimension: extreme categories get 1, middle categories get 0
  s_ers <- matrix(rep(0, categories), nrow = 1)
  s_ers[1, 1] <- 1  # Strongly disagree = extreme
  s_ers[1, categories] <- 1  # Strongly agree = extreme

  # Combine into full s-matrix
  s_matrix <- rbind(s_content, s_ers)

  # Generate m values for varying item difficulties
  # m values equally spaced between -0.5 and 0.5
  m_values <- seq(-0.5, 0.5, length.out = n_items)

  # Initialize storage matrices
  total_n <- N * n_groups
  obs_matrix <- matrix(NA, nrow = total_n, ncol = n_items)
  group_membership <- rep(NA, total_n)
  true_theta <- matrix(NA, nrow = total_n, ncol = 2)
  colnames(true_theta) <- c("theta", "ers")

  # Generate data for each group
  for (g in 1:n_groups) {

    cat("  Generating data for group", g, "...\n")

    # Row indices for this group
    group_rows <- ((g - 1) * N + 1):(g * N)

    # Draw latent traits from multivariate normal
    mu_theta <- c(theta_df$theta_mean[g], theta_df$ers_mean[g])
    sigma_theta <- diag(c(theta_df$theta_sd[g]^2, theta_df$ers_sd[g]^2))

    theta_matrix <- mvrnorm(n = N, mu = mu_theta, Sigma = sigma_theta)

    # Ensure no extreme values (truncate at ±3)
    theta_matrix[theta_matrix > 3] <- 3
    theta_matrix[theta_matrix < -3] <- -3

    # Store true theta values
    true_theta[group_rows, ] <- theta_matrix

    # Generate responses for each item
    for (j in 1:n_items) {

      # Item-specific thresholds (add m value for variation)
      item_thresholds <- thresholds + m_values[j]

      # Convert thresholds to intercepts
      # Following MNRM parameterization - generalized for any number of categories
      intercept <- rep(NA, categories)
      intercept[1] <- 0  # First category is reference

      # Build cumulative intercepts for remaining categories
      for (k in 2:categories) {
        if (k == 2) {
          intercept[k] <- -item_thresholds[k - 1] * alpha[1]
        } else {
          intercept[k] <- intercept[k - 1] - item_thresholds[k - 1] * alpha[1]
        }
      }

      # Calculate alpha-weighted s-matrix
      # For each dimension d and category k: alpha[d] * s_matrix[d, k]
      # Result is 2×categories matrix
      alpha_matrix <- diag(alpha) %*% s_matrix  # (2×2) %*% (2×categories) = (2×categories)

      # Calculate category probabilities for each person
      # P(category k) = exp(theta' * alpha_k + intercept_k) / sum over all categories

      # Initialize probability matrix
      prob_matrix <- matrix(NA, nrow = N, ncol = categories)

      for (person in 1:N) {

        # Linear predictor for each category
        # theta is 1×2, alpha_matrix is 2×4, result is 1×4
        linear_pred <- theta_matrix[person, ] %*% alpha_matrix + intercept

        # Convert to probabilities using softmax
        numerator <- exp(linear_pred)
        denominator <- sum(numerator)
        prob_matrix[person, ] <- numerator / denominator
      }

      # Sample responses based on probabilities
      obs_matrix[group_rows, j] <- apply(prob_matrix, 1, function(probs) {
        sample(1:categories, size = 1, prob = probs)
      })
    }

    # Store group membership
    group_membership[group_rows] <- g
  }

  # Create data frame
  item_names <- paste0("Item", 1:n_items)
  colnames(obs_matrix) <- item_names

  data <- as.data.frame(obs_matrix)
  data$country <- paste0("Group_", group_membership)
  data$true_theta <- true_theta[, 1]
  data$true_ers <- true_theta[, 2]

  # Return results
  result <- list(
    data = data,
    parameters = list(
      N = N,
      n_items = n_items,
      n_groups = n_groups,
      categories = categories,
      alpha = alpha,
      thresholds = thresholds,
      m_values = m_values,
      s_matrix = s_matrix,
      theta_df = theta_df
    ),
    true_theta = true_theta,
    group_membership = group_membership
  )

  return(result)
}


#' Generate All Simulation Conditions
#'
#' Creates datasets for all simulation conditions combining ERS levels
#' and sample sizes.
#'
#' @param n_items Number of items (10 or 20)
#' @param thresholds Base thresholds (average or difficult)
#' @param n_replications Number of replications per condition
#' @param seed Random seed for reproducibility
#' @return List of datasets for all conditions
generate_all_conditions <- function(n_items = 10,
                                    thresholds = c(-1, 0, 1),
                                    n_replications = 500,
                                    seed = 12345) {

  set.seed(seed)

  # Define ERS conditions
  ers_conditions <- list(
    low_ers = data.frame(
      theta_mean = c(0, 0),
      theta_sd = c(1, 1),
      ers_mean = c(0, 0),
      ers_sd = c(1, 1)
    ),
    medium_ers = data.frame(
      theta_mean = c(0, 0),
      theta_sd = c(1, 1),
      ers_mean = c(0, 0.5),
      ers_sd = c(1, 1)
    ),
    high_ers_diff = data.frame(
      theta_mean = c(0, 0),
      theta_sd = c(1, 1),
      ers_mean = c(0, 1.0),
      ers_sd = c(1, 1)
    ),
    equal_ers = data.frame(
      theta_mean = c(0, 0.3),
      theta_sd = c(1, 1),
      ers_mean = c(0.5, 0.5),
      ers_sd = c(1, 1)
    )
  )

  # Define sample sizes
  sample_sizes <- c(100, 500, 1000)

  # Generate all combinations
  all_datasets <- list()

  for (ers_name in names(ers_conditions)) {
    for (n_size in sample_sizes) {

      condition_name <- paste0(ers_name, "_N", n_size)
      cat("\n=== Generating condition:", condition_name, "===\n")

      # Generate multiple replications
      replications <- vector("list", n_replications)

      for (rep in 1:n_replications) {
        if (rep %% 50 == 0) {
          cat("  Replication", rep, "/", n_replications, "\n")
        }

        replications[[rep]] <- generate_mnrm_falk(
          N = n_size,
          n_items = n_items,
          alpha = c(1.5, 1.5),
          thresholds = thresholds,
          theta_df = ers_conditions[[ers_name]]
        )
      }

      all_datasets[[condition_name]] <- replications
    }
  }

  return(all_datasets)
}


#' Quick Test of Data Generation
#'
#' Generates a small test dataset to verify the function works
test_generation <- function() {

  cat("Testing MNRM-Falk data generation...\n")

  theta_df <- data.frame(
    theta_mean = c(0, 0),
    theta_sd = c(1, 1),
    ers_mean = c(0, 0.5),
    ers_sd = c(1, 1)
  )

  test_data <- generate_mnrm_falk(
    N = 100,
    n_items = 10,
    alpha = c(1.5, 1.5),
    thresholds = c(-1, 0, 1),
    theta_df = theta_df
  )

  cat("\nData structure:\n")
  print(str(test_data))

  cat("\nFirst 6 rows:\n")
  print(head(test_data$data))

  cat("\nGroup means (true theta):\n")
  print(aggregate(true_theta ~ country, data = test_data$data, mean))

  cat("\nGroup means (true ERS):\n")
  print(aggregate(true_ers ~ country, data = test_data$data, mean))

  cat("\nTest completed successfully!\n")

  return(test_data)
}

# Uncomment to test:
# test_result <- test_generation()
