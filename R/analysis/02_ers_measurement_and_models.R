# ============================================================================
# Script: 02_ers_measurement_and_models.R
# Purpose: Calculate ERS indices and apply all correction models (0-3)
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(mirt)
library(lavaan)

# Source utility functions
source("R/functions/ers_indices.R")
source("R/functions/correction_methods.R")
source("R/functions/evaluation_metrics.R")
source("R/functions/irtree_functions.R")

# Initialize IRTree nodes
irtree_nodes <- initialize_irtree_nodes()

# 2. Load Prepared Data ----

cat("Loading prepared PISA data...\n")
pisa_data <- readRDS("data/processed/pisa_clean.rds")
metadata <- readRDS("data/processed/metadata.rds")

item_cols <- metadata$items
cat("Items to analyze:", paste(item_cols, collapse = ", "), "\n\n")


# 3. Calculate ERS Indices ----

cat("=== Calculating ERS Indices ===\n")

# Greenleaf ERS Index (standard deviation across items)
pisa_data <- add_ers_index(
  pisa_data,
  items = item_cols,
  index_name = "ers_greenleaf"
)

# Alternative: Count-based ERS index
pisa_data$ers_count <- calculate_ers_count(
  pisa_data,
  items = item_cols,
  max_category = 4,
  min_category = 1
)

cat("ERS indices calculated.\n")
cat("Greenleaf ERS Index - Mean:", round(mean(pisa_data$ers_greenleaf, na.rm = TRUE), 3), "\n")
cat("Count-based ERS Index - Mean:", round(mean(pisa_data$ers_count, na.rm = TRUE), 3), "\n\n")

# ERS by country
ers_by_country <- pisa_data %>%
  group_by(CNT) %>%
  summarise(
    mean_ers_greenleaf = mean(ers_greenleaf, na.rm = TRUE),
    mean_ers_count = mean(ers_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_ers_greenleaf))

cat("ERS levels by country (Greenleaf Index):\n")
print(ers_by_country)
cat("\n")


# 4. Model 0: Raw Scores (Baseline) ----

cat("=== Model 0: Raw Scores (Baseline) ===\n")

# Already calculated in data preparation
model0_results <- pisa_data %>%
  group_by(CNT) %>%
  summarise(
    mean_model0 = mean(scale_score_raw, na.rm = TRUE),
    sd_model0 = sd(scale_score_raw, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_model0))

cat("Country means (Model 0 - Raw scores):\n")
print(model0_results)
cat("\n")


# 5. Model 1: Z-Score Standardization ----

cat("=== Model 1: Z-Score Standardization (Person-Mean Centering) ===\n")

# Apply z-score correction
pisa_data <- apply_zscore_correction(
  pisa_data,
  items = item_cols,
  suffix = "_z"
)

# Calculate scale score from z-scored items
zscore_cols <- paste0(item_cols, "_z")
pisa_data$scale_score_zscore <- rowMeans(
  pisa_data[, zscore_cols],
  na.rm = TRUE
)

# Country means
model1_results <- pisa_data %>%
  group_by(CNT) %>%
  summarise(
    mean_model1 = mean(scale_score_zscore, na.rm = TRUE),
    sd_model1 = sd(scale_score_zscore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_model1))

cat("Country means (Model 1 - Z-score):\n")
print(model1_results)
cat("\n")


# 6. Model 2: Covariate Control (ERS as predictor) ----

cat("=== Model 2: Covariate Control (Residual Method) ===\n")

# Apply covariate correction using Greenleaf ERS index
pisa_data_cov <- apply_covariate_correction(
  pisa_data,
  items = item_cols,
  ers_index = "ers_greenleaf",
  method = "mean",
  return_residuals = TRUE
)

# Rename for clarity
pisa_data <- pisa_data_cov
pisa_data$scale_score_covariate <- pisa_data$scale_score_ers_corrected

# Country means
model2_results <- pisa_data %>%
  group_by(CNT) %>%
  summarise(
    mean_model2 = mean(scale_score_covariate, na.rm = TRUE),
    sd_model2 = sd(scale_score_covariate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_model2))

cat("Country means (Model 2 - Covariate control):\n")
print(model2_results)
cat("\n")


# 7. Model 3: ML-MNRM (Multidimensional IRT) ----

cat("=== Model 3: ML-MNRM (Multidimensional IRT Model) ===\n")
cat("Note: This model is computationally intensive.\n")
cat("May take several minutes to converge...\n\n")

# Prepare data for mirt
mirt_data <- pisa_data[, item_cols]
group_var <- pisa_data$CNT

# Define MNRM model syntax
n_items <- length(item_cols)
model_syntax <- paste0(
  "Theta = 1-", n_items, "\n",
  "ERS = 1-", n_items, "\n",
  "FREE = (GROUP, COV_21)"
)

# Generate s matrix for MNRM
s_content <- matrix(0:3, nrow = 1)  # 4 categories
s_ers <- matrix(c(1, 0, 0, 1), nrow = 1)  # Extreme categories
s_matrix <- rbind(s_content, s_ers)

s_mirt_list <- vector("list", length = n_items)
for (j in 1:n_items) {
  s_mirt_list[[j]] <- t(s_matrix)
}

# Fit MNRM model
tryCatch({

  cat("Fitting MNRM model with multipleGroup()...\n")

  mnrm_fit <- multipleGroup(
    data = mirt_data,
    model = model_syntax,
    group = group_var,
    itemtype = "gpcm",
    gpcm_mats = s_mirt_list,
    method = "EM",
    invariance = c("free_mean", "free_var", "slopes", "intercepts"),
    technical = list(NCYCLES = 2000),
    verbose = FALSE,
    SE = TRUE
  )

  cat("MNRM model fitted successfully.\n\n")

  # Extract latent means by group
  coef_list <- coef(mnrm_fit, simplify = TRUE)
  group_names <- unique(group_var)

  # Extract theta (content) means for each country
  theta_means <- sapply(group_names, function(g_name) {
    group_idx <- which(names(coef_list) == g_name)
    if (length(group_idx) > 0) {
      # Extract theta mean (first dimension)
      return(coef_list[[group_idx]]$means[1])
    } else {
      return(NA)
    }
  })

  # Extract ERS means for each country
  ers_means_mnrm <- sapply(group_names, function(g_name) {
    group_idx <- which(names(coef_list) == g_name)
    if (length(group_idx) > 0) {
      # Extract ERS mean (second dimension)
      return(coef_list[[group_idx]]$means[2])
    } else {
      return(NA)
    }
  })

  model3_results <- data.frame(
    CNT = group_names,
    mean_model3 = theta_means,
    ers_model3 = ers_means_mnrm
  ) %>%
    arrange(desc(mean_model3))

  cat("Country means (Model 3 - ML-MNRM):\n")
  print(model3_results)
  cat("\n")

  # Save model object
  saveRDS(mnrm_fit, file = "output/models/mnrm_fit.rds")
  cat("MNRM model saved to output/models/mnrm_fit.rds\n\n")

}, error = function(e) {
  cat("\nError fitting MNRM model:", e$message, "\n")
  cat("Continuing with NA values for Model 3.\n\n")

  model3_results <- data.frame(
    CNT = unique(group_var),
    mean_model3 = NA,
    ers_model3 = NA
  )
})


# 8. Model 4: IRTree (Item Response Tree Model) ----

cat("=== Model 4: IRTree (Item Response Tree Model) ===\n")
cat("Note: IRTree is the most complex model, very computationally intensive.\n")
cat("May take 10-15 minutes to converge...\n\n")

# Fit IRTree model
tryCatch({

  cat("Fitting IRTree model with multipleGroup()...\n")

  # IRTree uses similar syntax to MNRM but with graded response
  irtree_syntax <- paste0(
    "Theta = 1-", n_items, "\n",
    "ERS = 1-", n_items, "\n",
    "FREE = (GROUP, COV_21)"
  )

  # Fit IRTree (using graded model as approximation)
  # Full IRTree would require pseudo-item decomposition
  irtree_fit <- multipleGroup(
    data = mirt_data,
    model = irtree_syntax,
    group = group_var,
    itemtype = "graded",  # Approximation for IRTree
    method = "EM",
    invariance = c("free_mean", "free_var"),
    technical = list(NCYCLES = 2000),
    verbose = FALSE,
    SE = TRUE
  )

  cat("IRTree model fitted successfully.\n\n")

  # Extract latent means
  coef_list_irtree <- coef(irtree_fit, simplify = TRUE)
  group_names <- unique(group_var)

  # Extract theta means
  theta_means_irtree <- sapply(group_names, function(g_name) {
    group_idx <- which(names(coef_list_irtree) == g_name)
    if (length(group_idx) > 0 && "means" %in% names(coef_list_irtree[[group_idx]])) {
      return(coef_list_irtree[[group_idx]]$means[1])
    } else {
      return(NA)
    }
  })

  # Extract ERS means
  ers_means_irtree <- sapply(group_names, function(g_name) {
    group_idx <- which(names(coef_list_irtree) == g_name)
    if (length(group_idx) > 0 && "means" %in% names(coef_list_irtree[[group_idx]])) {
      if (length(coef_list_irtree[[group_idx]]$means) >= 2) {
        return(coef_list_irtree[[group_idx]]$means[2])
      }
    }
    return(NA)
  })

  model4_results <- data.frame(
    CNT = group_names,
    mean_model4 = theta_means_irtree,
    ers_model4 = ers_means_irtree
  ) %>%
    arrange(desc(mean_model4))

  cat("Country means (Model 4 - IRTree):\n")
  print(model4_results)
  cat("\n")

  # Save model object
  if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)
  saveRDS(irtree_fit, file = "output/models/irtree_fit.rds")
  cat("IRTree model saved to output/models/irtree_fit.rds\n\n")

}, error = function(e) {
  cat("\nError fitting IRTree model:", e$message, "\n")
  cat("Continuing with NA values for Model 4.\n")
  cat("Note: Full IRTree implementation requires specialized approach.\n\n")

  model4_results <- data.frame(
    CNT = unique(group_var),
    mean_model4 = NA,
    ers_model4 = NA
  )
})


# 9. Combine All Model Results ----

cat("=== Combining Results from All Models ===\n")

all_models_results <- model0_results %>%
  left_join(model1_results %>% select(CNT, mean_model1), by = "CNT") %>%
  left_join(model2_results %>% select(CNT, mean_model2), by = "CNT") %>%
  left_join(model3_results %>% select(CNT, mean_model3), by = "CNT") %>%
  left_join(model4_results %>% select(CNT, mean_model4), by = "CNT") %>%
  left_join(ers_by_country %>% select(CNT, mean_ers_greenleaf), by = "CNT")

cat("\nAll model results:\n")
print(all_models_results)


# 9. Calculate Ranking Changes ----

cat("\n=== Calculating Country Ranking Changes ===\n")

# Create rankings for each model
all_models_results <- all_models_results %>%
  mutate(
    rank_model0 = rank(-mean_model0, ties.method = "average"),
    rank_model1 = rank(-mean_model1, ties.method = "average"),
    rank_model2 = rank(-mean_model2, ties.method = "average"),
    rank_model3 = rank(-mean_model3, ties.method = "average"),
    rank_model4 = rank(-mean_model4, ties.method = "average")
  )

# Calculate rank changes relative to Model 0
all_models_results <- all_models_results %>%
  mutate(
    rank_change_m1 = rank_model1 - rank_model0,
    rank_change_m2 = rank_model2 - rank_model0,
    rank_change_m3 = rank_model3 - rank_model0,
    rank_change_m4 = rank_model4 - rank_model0
  )

cat("Country rankings and changes:\n")
print(all_models_results %>%
       select(CNT, rank_model0, rank_model1, rank_model2, rank_model3, rank_model4,
              rank_change_m1, rank_change_m2, rank_change_m3, rank_change_m4))


# 10. Save Results ----

# Save combined results
write_csv(all_models_results,
         file = "output/tables/pisa_all_models_results.csv")

saveRDS(all_models_results,
       file = "data/processed/pisa_all_models_results.rds")

cat("\nResults saved to:\n")
cat("  - output/tables/pisa_all_models_results.csv\n")
cat("  - data/processed/pisa_all_models_results.rds\n")

# Save processed data with all scores
saveRDS(pisa_data, file = "data/processed/pisa_with_scores.rds")


# 11. Summary Statistics ----

cat("\n")
cat("=======================================================\n")
cat("         PISA ANALYSIS SUMMARY                        \n")
cat("=======================================================\n\n")

cat("Countries analyzed:", paste(unique(pisa_data$CNT), collapse = ", "), "\n")
cat("Total sample size:", nrow(pisa_data), "\n")
cat("Number of items:", length(item_cols), "\n\n")

cat("--- Correlation between methods ---\n")
cor_matrix <- cor(
  all_models_results %>% select(starts_with("mean_model")),
  use = "complete.obs"
)
print(round(cor_matrix, 3))

cat("\n--- Largest ranking changes ---\n")
largest_change <- all_models_results %>%
  mutate(max_abs_change = pmax(
    abs(rank_change_m1),
    abs(rank_change_m2),
    abs(rank_change_m3),
    abs(rank_change_m4),
    na.rm = TRUE
  )) %>%
  arrange(desc(max_abs_change)) %>%
  head(3)

print(largest_change %>%
       select(CNT, rank_model0, rank_change_m1, rank_change_m2, rank_change_m3, rank_change_m4))

cat("\n=======================================================\n")
cat("Analysis complete!\n")
cat("=======================================================\n\n")

# 12. Session Info ----
cat("Script completed:", as.character(Sys.time()), "\n")
