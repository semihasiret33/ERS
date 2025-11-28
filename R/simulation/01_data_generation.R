# ============================================================================
# Script: 01_data_generation.R
# Purpose: Generate simulated data using MNRM for different ERS conditions
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(mirt)

set.seed(12345)  # For reproducibility

# 2. Parameters ----

# Simulation conditions
N_REPLICATIONS <- 1000  # Number of simulation replications
N_COUNTRIES <- 2        # Number of countries/groups
N_PER_COUNTRY <- 500    # Sample size per country
N_ITEMS <- 10           # Number of Likert items
N_CATEGORIES <- 4       # Number of response categories (1-4)

# Item parameters
ITEM_ALPHA <- c(1, 1)        # Discrimination parameters for content and ERS
ITEM_THRESHOLDS <- c(-1, 0, 1)  # Threshold parameters

# ERS conditions to simulate
ers_conditions <- list(
  low_ers = data.frame(
    country = c("A", "B"),
    theta_mean = c(0, 0),
    theta_sd = c(1, 1),
    ers_mean = c(0, 0),      # Both countries low ERS
    ers_sd = c(0.5, 0.5)
  ),

  medium_ers = data.frame(
    country = c("A", "B"),
    theta_mean = c(0, 0),
    theta_sd = c(1, 1),
    ers_mean = c(0, 0.5),    # Country B medium ERS
    ers_sd = c(0.5, 0.75)
  ),

  high_ers_diff = data.frame(
    country = c("A", "B"),
    theta_mean = c(0, 0),      # Same true theta
    theta_sd = c(1, 1),
    ers_mean = c(0, 1.0),      # Country B high ERS
    ers_sd = c(0.5, 1.0)
  ),

  equal_ers = data.frame(
    country = c("A", "B"),
    theta_mean = c(0, 0.3),    # Different true theta
    theta_sd = c(1, 1),
    ers_mean = c(0.5, 0.5),    # Equal ERS
    ers_sd = c(0.75, 0.75)
  )
)


# 3. Functions ----

#' Generate MNRM Data for Multiple Countries
#'
#' Generates Likert response data under the MNRM (Falk & Cai, 2016) model
#' with separate content and ERS dimensions
#'
#' @param theta_df Data frame with columns: country, theta_mean, theta_sd,
#'   ers_mean, ers_sd
#' @param n_per_country Sample size per country
#' @param n_items Number of items
#' @param alpha Item discrimination parameters (length 2: content, ERS)
#' @param thresholds Item threshold parameters
#' @param categories Number of response categories
#' @return List containing: data (response matrix), true_theta (latent traits),
#'   country (group membership), true_country_means
generate_mnrm_data <- function(theta_df, n_per_country, n_items,
                              alpha = c(1, 1),
                              thresholds = c(-1, 0, 1),
                              categories = 4) {

  n_countries <- nrow(theta_df)
  total_n <- n_per_country * n_countries

  # Generate s matrix (scoring matrix)
  # Content dimension: 0, 1, 2, 3 (linear)
  s_content <- matrix(0:(categories - 1), nrow = 1)

  # ERS dimension: 1, 0, 0, 1 (extreme categories)
  s_ers <- matrix(c(1, 0, 0, 1), nrow = 1)

  # Combine into single scoring matrix
  s_matrix <- rbind(s_content, s_ers)

  # Initialize storage
  response_data <- matrix(NA, nrow = total_n, ncol = n_items)
  true_theta <- matrix(NA, nrow = total_n, ncol = 2)
  country_membership <- rep(NA, total_n)

  # Item difficulties (spread across range)
  item_difficulties <- seq(-0.5, 0.5, length = n_items)

  # Generate data for each country
  for (g in 1:n_countries) {

    # Indices for this country
    idx_start <- 1 + (g - 1) * n_per_country
    idx_end <- g * n_per_country
    idx <- idx_start:idx_end

    # Generate latent traits
    theta_content <- rnorm(n_per_country,
                          mean = theta_df$theta_mean[g],
                          sd = theta_df$theta_sd[g])
    theta_ers <- rnorm(n_per_country,
                      mean = theta_df$ers_mean[g],
                      sd = theta_df$ers_sd[g])

    # Store true theta
    true_theta[idx, ] <- cbind(theta_content, theta_ers)

    # Generate responses for each item
    for (j in 1:n_items) {

      # Calculate item intercepts from thresholds
      item_threshold <- thresholds - item_difficulties[j]
      intercepts <- c(
        0,
        -item_threshold[1] * alpha[1],
        -item_threshold[1] * alpha[1] - item_threshold[2] * alpha[1],
        -item_threshold[1] * alpha[1] - item_threshold[2] * alpha[1] -
          item_threshold[3] * alpha[1]
      )

      # Calculate slope parameters
      slopes <- alpha * s_matrix  # 2 x 4 matrix

      # Calculate category probabilities for each person
      # theta is n_per_country x 2, slopes is 2 x 4
      linear_pred <- cbind(theta_content, theta_ers) %*% slopes +
        matrix(intercepts, nrow = n_per_country, ncol = categories, byrow = TRUE)

      # Apply softmax
      exp_pred <- exp(linear_pred)
      cat_probs <- exp_pred / rowSums(exp_pred)

      # Sample responses
      responses <- apply(cat_probs, 1, function(probs) {
        sample(1:categories, size = 1, prob = probs)
      })

      response_data[idx, j] <- responses
    }

    # Store country membership
    country_membership[idx] <- theta_df$country[g]
  }

  # Add column names
  colnames(response_data) <- paste0("Item", 1:n_items)
  colnames(true_theta) <- c("true_theta", "true_ers")

  # Calculate true country means
  true_country_means <- theta_df %>%
    select(country, theta_mean, ers_mean)

  # Combine into data frame
  sim_data <- data.frame(
    country = country_membership,
    true_theta[, "true_theta"],
    true_theta[, "true_ers"],
    response_data
  )
  colnames(sim_data)[2:3] <- c("true_theta", "true_ers")

  result <- list(
    data = sim_data,
    true_country_means = true_country_means,
    parameters = list(
      alpha = alpha,
      thresholds = thresholds,
      item_difficulties = item_difficulties,
      s_matrix = s_matrix
    )
  )

  return(result)
}


# 4. Run Simulations ----

# Create directory to save results
if (!dir.exists("data/simulated")) {
  dir.create("data/simulated", recursive = TRUE)
}

# Generate one dataset per condition (for testing/development)
cat("Generating simulated datasets for each ERS condition...\n")

simulation_datasets <- list()

for (condition_name in names(ers_conditions)) {

  cat("  Generating:", condition_name, "\n")

  # Generate data
  sim_result <- generate_mnrm_data(
    theta_df = ers_conditions[[condition_name]],
    n_per_country = N_PER_COUNTRY,
    n_items = N_ITEMS,
    alpha = ITEM_ALPHA,
    thresholds = ITEM_THRESHOLDS,
    categories = N_CATEGORIES
  )

  # Store
  simulation_datasets[[condition_name]] <- sim_result

  # Save individual condition
  saveRDS(sim_result,
         file = paste0("data/simulated/", condition_name, "_dataset.rds"))
}

# Save all conditions together
saveRDS(simulation_datasets,
       file = "data/simulated/all_conditions.rds")

cat("\nSimulated datasets saved to data/simulated/\n")


# 5. Preview Results ----

cat("\n=== Preview of Simulated Data ===\n")
cat("\nCondition: high_ers_diff\n")
preview <- simulation_datasets$high_ers_diff

cat("\nTrue country parameters:\n")
print(preview$true_country_means)

cat("\nFirst few rows of simulated data:\n")
print(head(preview$data))

cat("\nObserved country means (raw item averages):\n")
observed_means <- preview$data %>%
  group_by(country) %>%
  summarise(
    mean_theta_true = mean(true_theta),
    mean_ers_true = mean(true_ers),
    mean_item_score = mean(c_across(starts_with("Item"))),
    .groups = "drop"
  )
print(observed_means)

cat("\nNote: mean_item_score is affected by ERS and will differ from true_theta\n")
cat("This is what we aim to correct with ERS adjustment methods.\n")


# 6. Session Info ----
cat("\n=== Session Info ===\n")
cat("Script completed:", as.character(Sys.time()), "\n")
